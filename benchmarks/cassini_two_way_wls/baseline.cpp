#include "../cassini_trajectory/CassiniForceModel.hpp"
#include "CsvTable.hpp"
#include "RKF45Integrator.hpp"
#include "DP853Integrator.hpp"
#include "dynamics/EphemerisInterpolator.hpp"
#include "dynamics/SpiceKernelSet.hpp"
#include "dynamics/SpiceStateProvider.hpp"
#include "filters/BatchLeastSquaresDriver.hpp"
#include "observations/radiometric/TwoWayMeasurement.hpp"
#include "observations/radiometric/VlbiMeasurement.hpp"
#include "observations/radiometric/PassRangeWhitening.hpp"
#include "observations/radiometric/SpiceLightTime.hpp"
#include "observations/radiometric/PropagationCorrections.hpp"
#include "simulation/CassiniRadioKernels.hpp"
#include <chrono>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <numbers>
#include <random>
#include <array>
#include <sstream>

namespace radio=fd::observations::radiometric;
using State=Eigen::Matrix<double,6,1>;
using Parameters=Eigen::Matrix<double,8,1>;
using Matrix=Eigen::Matrix<double,8,8>;
using Sensitivity=Eigen::Matrix<double,6,8>;
namespace {
struct Options {
    std::string utc="2004-10-09T00:00:00",truth="simulation";
    std::filesystem::path output="results/cassini_two_way_wls";
    std::filesystem::path observationInput;
    double days=14,cadence=600,count=60,mass=3220,area=12,cr=1.42,rangeSigma=.01,dopplerSigma=.01,maxStep=60;
    Eigen::Vector3d recoil{-2.725e-13,-2.725e-13,-5.232237415805135e-12};
    int vlbiCount=0;double vlbiSigma=1e-5,initialGuessScale=1;
    double rangeBiasSigma=.01,rangeCadence=1200,passHours=6;
    int iterations=10;unsigned seed=2026;bool noise=true,passSchedule=false,exportCovariance=false,generateOnly=false;
};
Options parse(int argc,char** argv) {
    Options o;
    for (int i=1;i<argc;++i) {
        std::string key=argv[i];
        if (key=="--generate-only") {o.generateOnly=true;continue;}
        if (key=="--export-covariance") {o.exportCovariance=true;continue;}
        if (key=="--pass-schedule") {o.passSchedule=true;o.cadence=60;continue;}
        if (key=="--no-noise") {o.noise=false;continue;}
        if (key=="--help") {std::cout<<"--start-utc UTC --days D --cadence-s S --count-s S --mass-kg M --area-m2 A --cr CR --thermal-x-km-s2 X --thermal-y-km-s2 Y --thermal-z-km-s2 Z --output DIR --truth simulation --no-noise --seed N --max-iterations N --max-step-s S --range-sigma-m M --doppler-sigma-hz H --vlbi-count N --vlbi-sigma-m M --initial-guess-scale F --pass-schedule --range-bias-sigma-m M --range-cadence-s S --pass-hours H --export-covariance --generate-only --observation-input DIR\n";std::exit(0);}
        if (++i==argc) throw std::invalid_argument("Missing value for "+key);
        if (key=="--output") o.output=argv[i];
        else if (key=="--observation-input") o.observationInput=argv[i];
        else if (key=="--start-utc") o.utc=argv[i];
        else if (key=="--truth") o.truth=argv[i];
        else {
            const double value=std::stod(argv[i]);if (!std::isfinite(value)) throw std::invalid_argument("Nonfinite option");
            if (key=="--days") o.days=value;else if (key=="--cadence-s") o.cadence=value;
            else if (key=="--count-s") o.count=value;else if (key=="--mass-kg") o.mass=value;
            else if (key=="--area-m2") o.area=value;else if (key=="--cr") o.cr=value;
            else if (key=="--thermal-x-km-s2") o.recoil[0]=value;else if (key=="--thermal-y-km-s2") o.recoil[1]=value;
            else if (key=="--thermal-z-km-s2") o.recoil[2]=value;else if (key=="--range-sigma-m") o.rangeSigma=value/1000;
            else if (key=="--doppler-sigma-hz") o.dopplerSigma=value;else if (key=="--max-step-s") o.maxStep=value;
            else if (key=="--vlbi-count" && value>=0 && value<=1000 && value==std::floor(value)) o.vlbiCount=static_cast<int>(value);
            else if (key=="--vlbi-sigma-m" && value>0) o.vlbiSigma=value/1000;
            else if (key=="--initial-guess-scale" && value>0) o.initialGuessScale=value;
            else if (key=="--range-bias-sigma-m" && value>=0) o.rangeBiasSigma=value/1000;
            else if (key=="--range-cadence-s" && value>0) o.rangeCadence=value;
            else if (key=="--pass-hours" && value>0 && value<=8) o.passHours=value;
            else if (key=="--seed" && value>=0 && value==std::floor(value) && value<=4294967295.) o.seed=static_cast<unsigned>(value);
            else if (key=="--max-iterations" && value>=2 && value<=30 && value==std::floor(value)) o.iterations=static_cast<int>(value);
            else throw std::invalid_argument("Unknown/invalid option "+key);
        }
    }
    if (!(o.days>0 && o.days<=14 && o.cadence>0 && o.count>0 && o.count<=o.cadence && o.count<=300 && o.mass>0 && o.area>0 && o.cr>=1 && o.cr<=5./3 && o.rangeSigma>0 && o.dopplerSigma>0 && o.maxStep>0 && o.maxStep<=120) || o.truth!="simulation") throw std::invalid_argument("Invalid physical/schedule options");
    if (o.passSchedule && (o.cadence!=60 || o.count!=60 || o.days!=std::floor(o.days) || o.rangeCadence<o.cadence || std::fmod(o.rangeCadence,o.cadence)!=0)) throw std::invalid_argument("Pass mode requires whole days, 60 s counts/cadence and an integral range-cadence multiple");
    return o;
}
double epoch(const std::string& utc) {double et;str2et_c(utc.c_str(),&et);od::throwIfSpiceFailed("Epoch");return et;}
Parameters scales() {return (Parameters()<<10000,10000,10000,10,10,10,1,1).finished();}
class Arc final : public fd::dynamics::StateProvider {
    double t0_;const fd::dynamics::StateProvider& origin_;
    od::EphemerisInterpolator history_;
public:
    Parameters initial;
    Arc(double t0,double duration,const Parameters& x,const fd::simulation::CassiniForceModel& force,
        const fd::dynamics::StateProvider& origin,double maxStep):t0_(t0),origin_(origin),initial(x) {
        const Parameters scale=scales(),inv=scale.cwiseInverse();Eigen::VectorXd augmented=Eigen::VectorXd::Zero(72);
        augmented.head<8>()=inv.asDiagonal()*x;Eigen::Map<Matrix>(augmented.data()+8).setIdentity();
        od::RKF45Integrator::Options options;options.absoluteTolerance=1e-13;options.relativeTolerance=1e-13;options.initialStep=10;options.maximumStep=maxStep;
        const auto rate=[&](double elapsed,od::RKF45Integrator::ConstStateRef a,od::RKF45Integrator::StateRef da) {
            const Parameters state=scale.asDiagonal()*a.head<8>();const auto f=force.evaluate(t0+elapsed,state.head<3>(),state.segment<3>(3),true,state[6],state[7]);
            Parameters derivative=Parameters::Zero();derivative.head<6>()<<state.segment<3>(3),f.acceleration;da.head<8>()=inv.asDiagonal()*derivative;
            Matrix jac=Matrix::Zero();jac.block<3,3>(0,3).setIdentity();jac.block<3,3>(3,0)=f.positionJacobian;jac.block<3,3>(3,3)=f.velocityJacobian;
            jac.block<3,1>(3,6)=f.srpScaleJacobian;jac.block<3,1>(3,7)=f.thermalScaleJacobian;
            Eigen::Map<Matrix>(da.data()+8)=inv.asDiagonal()*jac*scale.asDiagonal()*Eigen::Map<const Matrix>(a.data()+8);
        };
        const auto result=od::RKF45Integrator(options).integrate(0,augmented,duration,rate);
        for (const auto& node:result.history) history_.addNode(node.tdb,node.state,node.derivative);
    }
    State relative(double t) const {return (scales().asDiagonal()*history_.interpolate(t-t0_).head<8>()).head<6>();}
    Sensitivity transition(double t) const {const auto x=history_.interpolate(t-t0_);const Matrix phi=scales().asDiagonal()*Eigen::Map<const Matrix>(x.data()+8)*scales().cwiseInverse().asDiagonal();return phi.topRows<6>();}
    fd::dynamics::CartesianState stateAt(fd::opnav::TdbEpoch t) const override {
        const auto x=relative(t.secondsPastJ2000);const auto saturn=origin_.stateAt(t);
        return {x.head<3>()+saturn.positionKm,x.tail<3>()+saturn.velocityKmPerSec};
    }
    std::string_view referenceFrame() const noexcept override {return "J2000";}
    std::string_view referenceOrigin() const noexcept override {return "SSB";}
};
class AffineArc final:public fd::dynamics::StateProvider {
    const Arc& arc_;Parameters correction_;
public:
    AffineArc(const Arc& arc,const Parameters& correction):arc_(arc),correction_(correction) {}
    fd::dynamics::CartesianState stateAt(fd::opnav::TdbEpoch t) const override {
        auto s=arc_.stateAt(t);const State delta=arc_.transition(t.secondsPastJ2000)*correction_;
        s.positionKm+=delta.head<3>();s.velocityKmPerSec+=delta.tail<3>();return s;
    }
    std::string_view referenceFrame() const noexcept override {return "J2000";}
    std::string_view referenceOrigin() const noexcept override {return "SSB";}
};
class SavedSimulation final:public fd::dynamics::StateProvider {
    od::EphemerisInterpolator history_;const fd::dynamics::StateProvider& origin_;
public:
    SavedSimulation(const std::filesystem::path& path,const fd::simulation::CassiniForceModel& force,const fd::dynamics::StateProvider& origin):origin_(origin) {
        std::ifstream input(path);if (!input) throw std::runtime_error("Cannot read retained simulation CSV");
        std::string line;std::getline(input,line);
        if (line.find("epoch_tdb,utc,model,x_km,y_km,z_km,vx_km_s,vy_km_s,vz_km_s")==std::string::npos) throw std::runtime_error("Unexpected retained trajectory CSV schema");
        while (std::getline(input,line)) {
            std::stringstream row(line);std::vector<std::string> fields;std::string field;while (std::getline(row,field,',')) fields.push_back(field);
            if (fields.size()<10) throw std::runtime_error("Malformed retained trajectory row");
            if (fields[3]!="multibody_j2_srp_thermal") continue;
            const double t=std::stod(fields[1]);State x,dx;for (int j=0;j<6;++j) x[j]=std::stod(fields[4+j]);
            dx<<x.tail<3>(),force.evaluate(t,x.head<3>(),x.tail<3>()).acceleration;history_.addNode(t,x,dx);
        }
        if (history_.history().size()<2) throw std::runtime_error("Missing saved dynamical trajectory");
    }
    State relative(double t) const {return history_.interpolate(t);}
    fd::dynamics::CartesianState stateAt(fd::opnav::TdbEpoch t) const override {
        const auto x=relative(t.secondsPastJ2000);const auto origin=origin_.stateAt(t);return {x.head<3>()+origin.positionKm,x.tail<3>()+origin.velocityKmPerSec};
    }
    std::string_view referenceFrame() const noexcept override {return "J2000";}
    std::string_view referenceOrigin() const noexcept override {return "SSB";}
};
struct VlbiObservation {double elapsed,elevationOne,elevationTwo,truth,noise,value;std::size_t pair;};
struct Observation {double elapsed,elevation;std::size_t station;Eigen::Vector2d truth,noise,value;int pass=-1;bool rangeActive=true;double rangeBias=0,rangeIndependentNoise=0;};
Eigen::Vector2d value(const radio::TwoWayLink& link,double t,double count,const fd::dynamics::StateProvider& station,const fd::dynamics::StateProvider& probe,const radio::FrequencyRamp& ramp,double referenceHz,const radio::SpiceTdbClock& clock) {
    const auto c=link.count(t,count,station,probe,ramp,referenceHz,clock);return {c.end.rangeKm(),c.residualHz};
}
std::pair<double,double> elevations(const std::string& name,const radio::TwoWaySolution& s) {
    return {radio::spiceStationElevationRadians(name,s.uplink.emissionTdb,s.uplink.receiver.positionKm-s.uplink.emitter.positionKm),radio::spiceStationElevationRadians(name,s.downlink.receptionTdb,s.downlink.emitter.positionKm-s.downlink.receiver.positionKm)};
}
template<class Derived> void jsonVector(std::ostream& stream,const Eigen::MatrixBase<Derived>& x) {stream<<'[';for (Eigen::Index i=0;i<x.size();++i) stream<<(i?",":"")<<x.derived().coeff(i);stream<<']';}
void checkForcesAndSTM(const fd::simulation::CassiniForceModel& force,const Arc& nominal,double t0,double duration) {
    for (double fraction:{.1,.5,.9}) {
        const double t=t0+fraction*duration;const auto x=nominal.relative(t);const auto f=force.evaluate(t,x.head<3>(),x.tail<3>(),true);Eigen::Matrix3d numerical;
        for (int j=0;j<3;++j) {State a=x,b=x;a[j]+=.5;b[j]-=.5;numerical.col(j)=force.evaluate(t,a.head<3>(),a.tail<3>()).acceleration-force.evaluate(t,b.head<3>(),b.tail<3>()).acceleration;}
        if ((numerical-f.positionJacobian).norm()/numerical.norm()>1e-6) throw std::runtime_error("Force partials disagree with finite differences");
    }
    const double dt=std::min(duration,86400.);Parameters perturbation;perturbation<<.1,-.07,.04,1e-7,-8e-8,6e-8,1,-1;
    od::DP853Integrator::Options options;options.absoluteTolerance=1e-13;options.relativeTolerance=1e-14;options.maximumStep=30;double step=10;
    const auto dynamics=[&](double elapsed,State x)->State {State dx;dx<<x.tail<3>(),force.evaluate(t0+elapsed,x.head<3>(),x.tail<3>(),false,nominal.initial[6]+perturbation[6],nominal.initial[7]+perturbation[7]).acceleration;return dx;};
    const State integrated=od::DP853Integrator(options).integrate(0,State((nominal.initial+perturbation).head<6>()),dt,dynamics,step);
    const State linear=nominal.relative(t0+dt)+nominal.transition(t0+dt)*perturbation;
    if ((integrated.head<3>()-linear.head<3>()).norm()*1000>.05) throw std::runtime_error("STM nonlinear propagation check failed");
    std::cout<<"Force Jacobians and one-day STM nonlinear check passed.\n"<<std::flush;
}
} // namespace
int main(int argc,char** argv) {
 try {
    const auto o=parse(argc,argv);const auto wall=std::chrono::steady_clock::now();const std::filesystem::path root(DEEPNAV_SOURCE_DIR);
    od::SpiceErrorModeGuard errors("RETURN","NONE");auto paths=fd::simulation::cassiniRadioKernels(root/"Kernels");
    for (const char* file:{"Cassini/fk/cas_v40.tf","Cassini/sclk/cas00172.tsc","Cassini/ck/04256_04292py_as_flown.bc","Cassini/ck/04292_04320py_as_flown.bc","Cassini/ck/04281_04286ra.bc","Cassini/ck/04286_04291ra.bc","Cassini/ck/04291_04296ra.bc","Cassini/ck/04296_04301ra.bc","Cassini/ck/04293_04298ra.bc"}) paths.push_back(root/"Kernels"/file);
    fd::dynamics::SpiceKernelSet kernels(paths);const CassiniAttitude attitude(paths);const Sat427Gravity gravity(root/"Kernels/Cassini/spk/200128RU_SCPSE_04251_04303.bsp.lbl");
    const double t0=epoch(o.utc),duration=o.days*86400;double radii[3];SpiceInt n;bodvrd_c("SATURN","RADII",3,&n,radii);od::throwIfSpiceFailed("Saturn radius");
    const fd::simulation::CassiniForceModel force(gravity,attitude,fd::simulation::cassiniPerturbers(gravity),2,o.mass,o.area,o.cr,o.recoil,epoch("2005-08-13T00:00:00"),radii[0]);
    fd::dynamics::SpiceStateProvider saturn("SATURN"),sun("SUN");const radio::SpiceTdbClock clock;
    const SavedSimulation truth(root/"results/cassini_two_way_wls/pass_schedule/truth/trajectory.csv",force,saturn);
    if (o.exportCovariance) {
        // Reuse the completed estimate; no observations or WLS results are changed.
        Parameters estimate;Matrix covariance;std::string line,field;
        std::ifstream input(o.output/"initial_estimate.csv");
        if (!input || !std::getline(input,line)) throw std::runtime_error("Missing saved estimate");
        for (int i=0;i<8;++i) {if (!std::getline(input,line)) throw std::runtime_error("Incomplete estimate");std::istringstream row(line);for (int j=0;j<5;++j) {if (!std::getline(row,field,',')) throw std::runtime_error("Invalid estimate row");}estimate[i]=std::stod(field);}
        std::ifstream cov(o.output/"initial_covariance.csv");
        for (int i=0;i<8;++i) {if (!std::getline(cov,line)) throw std::runtime_error("Incomplete covariance");std::istringstream row(line);for (int j=0;j<8;++j) {if (!std::getline(row,field,',')) throw std::runtime_error("Invalid covariance row");covariance(i,j)=std::stod(field);}}
        if (!estimate.allFinite() || !covariance.allFinite() || covariance.llt().info()!=Eigen::Success) throw std::runtime_error("Invalid saved solution");
        const Arc fitted(t0,duration,estimate,force,saturn,o.maxStep);
        std::ofstream out(o.output/"trajectory_uncertainty.csv");
        out<<std::setprecision(17)<<"elapsed_s,radial_sigma_m,along_sigma_m,cross_sigma_m,formal_position_sigma_m,position_error_m,formal_velocity_sigma_mm_s,velocity_error_mm_s\n";
        for (std::size_t k=0;;++k) {
            const double elapsed=std::min(k*o.cadence,duration),t=t0+elapsed;const State target=truth.relative(t);
            const Eigen::Vector3d radial=target.head<3>().normalized(),normal=target.head<3>().cross(target.tail<3>()).normalized(),along=normal.cross(radial);
            const Sensitivity phi=fitted.transition(t);const Eigen::Matrix<double,6,6> P=phi*covariance*phi.transpose();const Eigen::Matrix3d position=P.topLeftCorner<3,3>();
            out<<elapsed;for (const auto& direction:{radial,along,normal}) out<<','<<std::sqrt(std::max(0.,direction.dot(position*direction)))*1000;
            const State error=fitted.relative(t)-target;
            out<<','<<std::sqrt(position.trace())*1000<<','<<error.head<3>().norm()*1000<<','<<std::sqrt(P.bottomRightCorner<3,3>().trace())*1e6<<','<<error.tail<3>().norm()*1e6<<'\n';if (elapsed==duration) break;
        }
        std::cout<<"Exported propagated covariance from saved estimate.\n";return 0;
    }
    Parameters trueInitial;trueInitial<<truth.relative(t0),1,1;truth.relative(t0+duration);
    std::cout<<"Integrating retained nominal orbit and STM for "<<o.days<<" days.\n"<<std::flush;
    auto nominal=std::make_shared<Arc>(t0,duration,trueInitial,force,saturn,o.maxStep);checkForcesAndSTM(force,*nominal,t0,duration);
    const std::array<std::string,3> names{"DSS-14","DSS-43","DSS-63"};
    std::vector<std::unique_ptr<fd::dynamics::SpiceStateProvider>> stations;std::vector<radio::TwoWayLink> links;std::vector<radio::DelayModel> corrections;
    const radio::LinkConfig transponder{880./749,1e-6};const radio::TwoWayLink geometric(radio::libraryReceptionSolver(),transponder);
    const radio::FrequencyRamp ramp{t0,7.175e9,0};const double referenceHz=transponder.turnaroundRatio*ramp.frequencyHz;
    const double solarGm=fd::simulation::cassiniKernelGM("SUN");
    for (std::size_t j=0;j<names.size();++j) {
        stations.push_back(std::make_unique<fd::dynamics::SpiceStateProvider>(names[j]));const auto* station=stations.back().get();const auto name=names[j];
        const auto correction=[&,station,name](const radio::LegSolution& leg,const fd::dynamics::StateProvider& emitter,const auto&) {
            const bool up=&emitter==station;const Eigen::Vector3d direction=up?(leg.receiver.positionKm-leg.emitter.positionKm).eval():(leg.emitter.positionKm-leg.receiver.positionKm).eval();
            const double elevation=radio::spiceStationElevationRadians(name,up?leg.emissionTdb:leg.receptionTdb,direction);
            if (elevation<=0) throw std::runtime_error("Link below horizon");return radio::pointMassShapiroSeconds(leg,sun,solarGm)+radio::zenithPathDelaySeconds(2.3,elevation);
        };
        corrections.push_back(correction);links.emplace_back(radio::libraryReceptionSolver({},correction),transponder);
    }
    std::filesystem::create_directories(o.output);std::ofstream schedule(o.output/"schedule.csv");schedule<<std::setprecision(17)<<"elapsed_s,station,minimum_elevation_deg,tracking\n";
    std::mt19937_64 rng(o.seed);std::normal_distribution<double> gaussian;std::vector<Observation> observations;std::array<int,3> stationCounts{};int candidates=0;
    // Reception must follow the arc start sufficiently to cover both retarded
    // spacecraft epochs of the end-tagged count. No backward extrapolation.
    const double first=std::ceil(6000/o.cadence)*o.cadence;
    if (o.observationInput.empty() && !o.passSchedule) for (double elapsed=first;elapsed<=duration;elapsed+=o.cadence) {
        ++candidates;const double t=t0+elapsed;double best=-1e9;std::size_t chosen=0;
        for (std::size_t j=0;j<names.size();++j) {const auto a=elevations(names[j],geometric.solve(t-o.count,*stations[j],truth)),b=elevations(names[j],geometric.solve(t,*stations[j],truth));const double e=std::min({a.first,a.second,b.first,b.second});if (e>best) {best=e;chosen=j;}}
        const bool tracking=best>=10*std::numbers::pi/180;schedule<<elapsed<<','<<names[chosen]<<','<<best*180/std::numbers::pi<<','<<tracking<<'\n';
        if (!tracking) continue;
        const auto count=links[chosen].count(t,o.count,*stations[chosen],truth,ramp,referenceHz,clock);
        if (count.start.spacecraftReceptionTdb<t0 || count.end.spacecraftTransmissionTdb>t0+duration) throw std::runtime_error("Retarded count outside dynamical arc");
        Eigen::Vector2d noise=Eigen::Vector2d::Zero();if (o.noise) noise<<o.rangeSigma*gaussian(rng),o.dopplerSigma*gaussian(rng);
        Eigen::Vector2d y(count.end.rangeKm(),count.residualHz);observations.push_back({elapsed,best,chosen,y,noise,y+noise});++stationCounts[chosen];
    }
    std::ofstream passFile(o.output/"passes.csv");passFile<<std::setprecision(17)<<"pass_id,station,start_elapsed_s,end_elapsed_s,usable_hours,doppler_points,range_points,true_range_bias_m,minimum_elevation_deg\n";
    const int requestedSamples=static_cast<int>(std::round(o.passHours*3600/o.cadence)),rangeEvery=static_cast<int>(o.rangeCadence/o.cadence);
    std::mt19937_64 biasRng(static_cast<std::uint64_t>(o.seed)+2000003);std::normal_distribution<double> biasGaussian;
    double previousPassEnd=first-o.count;
    if (o.observationInput.empty() && o.passSchedule) for (int day=0;day<static_cast<int>(o.days);++day) {
        std::vector<double> times;std::array<std::vector<double>,3> elevation;
        for (double t=std::max(first,day*86400.-o.passHours*1800+o.count);t<=std::min(duration,(day+1)*86400.+o.passHours*1800);t+=o.cadence) {
            times.push_back(t);++candidates;
            for (std::size_t j=0;j<3;++j) {const auto a=elevations(names[j],geometric.solve(t0+t-o.count,*stations[j],truth)),b=elevations(names[j],geometric.solve(t0+t,*stations[j],truth));elevation[j].push_back(std::min({a.first,a.second,b.first,b.second}));}
        }
        std::size_t chosen=0,begin=0,length=0;double best=-1;
        for (std::size_t j=0;j<3;++j) {
            for (std::size_t a=0;a<times.size();) {
                if (elevation[j][a]<10*std::numbers::pi/180) {++a;continue;}
                std::size_t b=a;while (b<times.size() && elevation[j][b]>=10*std::numbers::pi/180) ++b;
                const std::size_t n=std::min(static_cast<std::size_t>(requestedSamples),b-a);
                double sum=0;for (std::size_t k=a;k<a+n;++k) sum+=elevation[j][k];
                for (std::size_t w=a;w+n<=b;++w) {
                    const double midpoint=(times[w]-o.count+times[w+n-1])/2;
                    if (midpoint<day*86400. || midpoint>=(day+1)*86400. || times[w]-o.count<previousPassEnd) {if (w+n<b) sum+=elevation[j][w+n]-elevation[j][w];continue;}
                    const double score=n*100.+(j==static_cast<std::size_t>(day%3)?100.*(requestedSamples+1):0.)+sum/n;
                    if (score>best) {best=score;chosen=j;begin=w;length=n;}
                    if (w+n<b) sum+=elevation[j][w+n]-elevation[j][w];
                }a=b;
            }
        }
        if (!length) throw std::runtime_error("No visible tracking pass for a day");
        const double bias=o.noise?o.rangeBiasSigma*biasGaussian(biasRng):0;int ranges=0;double minimum=2;
        for (std::size_t k=0;k<length;++k) {
            const double elapsed=times[begin+k],t=t0+elapsed;const auto count=links[chosen].count(t,o.count,*stations[chosen],truth,ramp,referenceHz,clock);
            if (count.start.spacecraftReceptionTdb<t0 || count.end.spacecraftTransmissionTdb>t0+duration) throw std::runtime_error("Pass count outside simulated arc");
            const bool hasRange=(k+1)%rangeEvery==0;Eigen::Vector2d noise=Eigen::Vector2d::Zero();
            if (o.noise) {noise[1]=o.dopplerSigma*gaussian(rng);if (hasRange) noise[0]=o.rangeSigma*gaussian(rng);}
            const double independent=noise[0];if (hasRange) {noise[0]+=bias;++ranges;}
            Eigen::Vector2d y(count.end.rangeKm(),count.residualHz);observations.push_back({elapsed,elevation[chosen][begin+k],chosen,y,noise,y+noise,day,hasRange,bias,independent});++stationCounts[chosen];minimum=std::min(minimum,elevation[chosen][begin+k]);
            schedule<<elapsed<<','<<names[chosen]<<','<<elevation[chosen][begin+k]*180/std::numbers::pi<<",1\n";
        }
        std::cout<<"Pass "<<day<<": "<<names[chosen]<<", "<<length*o.cadence/3600<<" usable hours, "<<ranges<<" ranges.\n"<<std::flush;
        previousPassEnd=times[begin+length-1];
        passFile<<day<<','<<names[chosen]<<','<<times[begin]-o.count<<','<<times[begin+length-1]<<','<<length*o.cadence/3600<<','<<length<<','<<ranges<<','<<bias*1000<<','<<minimum*180/std::numbers::pi<<'\n';
    }passFile.close();schedule.close();
    if (!o.observationInput.empty()) {
        if (!o.passSchedule || std::filesystem::absolute(o.observationInput).lexically_normal()==std::filesystem::absolute(o.output).lexically_normal()) throw std::invalid_argument("Observation import requires pass mode and a separate input directory");
        const cassini_benchmark::CsvTable table(o.observationInput/"observations.csv");double previous=0;
        for (const auto& row:table.rows) {
            Observation a;a.elapsed=table.number(row,"elapsed_s");a.elevation=10*std::numbers::pi/180;
            auto station=std::find(names.begin(),names.end(),row[table.column("station")]);if (station==names.end()) throw std::runtime_error("Unknown imported station");a.station=station-names.begin();
            a.pass=static_cast<int>(table.number(row,"pass_id"));a.rangeActive=table.number(row,"range_active")==1;
            if (!(a.elapsed>=first && a.elapsed<=duration && a.elapsed>previous) || a.pass<0 || a.pass>=static_cast<int>(o.days)) throw std::runtime_error("Invalid imported epoch/pass");previous=a.elapsed;
            a.truth<<table.number(row,"truth_range_km"),table.number(row,"truth_doppler_hz");
            a.noise<<(a.rangeActive?table.number(row,"range_noise_m")/1000:0),table.number(row,"doppler_noise_hz");
            a.value<<(a.rangeActive?table.number(row,"observed_range_km"):a.truth[0]),table.number(row,"observed_doppler_hz");
            a.rangeBias=table.number(row,"range_bias_m")/1000;a.rangeIndependentNoise=table.number(row,"range_independent_noise_m")/1000;
            if (!a.truth.allFinite() || !a.value.allFinite() || !a.noise.allFinite()) throw std::runtime_error("Nonfinite imported observation");
            observations.push_back(a);++stationCounts[a.station];++candidates;
        }
        for (const char* file:{"schedule.csv","passes.csv"}) std::filesystem::copy_file(o.observationInput/file,o.output/file,std::filesystem::copy_options::overwrite_existing);
    }
    if (observations.size()<10) throw std::runtime_error("Insufficient visible observations");
    const std::array<std::array<std::size_t,2>,3> pairs{{{0,1},{0,2},{1,2}}};
    std::vector<radio::VlbiLink> vlbiLinks;
    for (const auto pair:pairs) vlbiLinks.emplace_back(radio::libraryReceptionSolver({},corrections[pair[0]]),corrections[pair[1]]);
    const radio::VlbiLink vlbiGeometry(radio::libraryReceptionSolver());
    std::vector<VlbiObservation> eligible,vlbi;
    std::array<int,3> pairCounts{};
    if (o.observationInput.empty() && o.vlbiCount) {
        for (double elapsed=first;elapsed<=duration;elapsed+=o.cadence) for (std::size_t j=0;j<pairs.size();++j) {
            const auto pair=pairs[j];const auto g=vlbiGeometry.solve(t0+elapsed,truth,*stations[pair[0]],*stations[pair[1]]);
            const double a=radio::spiceStationElevationRadians(names[pair[0]],g.first.receptionTdb,g.first.emitter.positionKm-g.first.receiver.positionKm);
            const double b=radio::spiceStationElevationRadians(names[pair[1]],g.second.receptionTdb,g.second.emitter.positionKm-g.second.receiver.positionKm);
            if (std::min(a,b)>=10*std::numbers::pi/180) eligible.push_back({elapsed,a,b,0,0,0,j});
        }
        std::mt19937_64 vlbiRng(static_cast<std::uint64_t>(o.seed)+1000003);std::normal_distribution<double> vlbiGaussian;
        for (int k=0;k<o.vlbiCount;++k) {
            const double lo=first+(duration-first)*k/o.vlbiCount,hi=first+(duration-first)*(k+1)/o.vlbiCount,mid=(lo+hi)/2;
            auto best=eligible.end();double score=1e100;
            for (auto it=eligible.begin();it!=eligible.end();++it) {
                bool used=false;for (const auto& selected:vlbi) if (selected.elapsed==it->elapsed) used=true;
                if (used) continue;
                const double cost=std::abs(it->elapsed-mid)+(it->pair==static_cast<std::size_t>(k%3)?0:3600);
                if (cost<score) {score=cost;best=it;}
            }
            if (best==eligible.end()) throw std::runtime_error("Insufficient distinct jointly visible VLBI epochs");
            auto a=*best;const auto pair=pairs[a.pair];a.truth=vlbiLinks[a.pair].solve(t0+a.elapsed,truth,*stations[pair[0]],*stations[pair[1]]).differentialRangeKm;
            a.noise=o.noise?o.vlbiSigma*vlbiGaussian(vlbiRng):0;a.value=a.truth+a.noise;vlbi.push_back(a);++pairCounts[a.pair];
        }
    }
    if (!o.observationInput.empty() && o.vlbiCount) {
        const cassini_benchmark::CsvTable table(o.observationInput/"vlbi_observations.csv");
        if (table.rows.size()!=static_cast<std::size_t>(o.vlbiCount)) throw std::runtime_error("Imported VLBI count mismatch");
        for (const auto& row:table.rows) {
            std::size_t pair=0;for (;pair<pairs.size();++pair) if (row[table.column("station_one")]==names[pairs[pair][0]] && row[table.column("station_two")]==names[pairs[pair][1]]) break;
            if (pair==pairs.size()) throw std::runtime_error("Unknown imported VLBI pair");
            VlbiObservation a{table.number(row,"elapsed_s"),table.number(row,"elevation_one_deg")*std::numbers::pi/180,table.number(row,"elevation_two_deg")*std::numbers::pi/180,table.number(row,"truth_delay_km"),table.number(row,"noise_m")/1000,table.number(row,"observed_delay_km"),pair};
            if (!(a.elapsed>=first && a.elapsed<=duration && a.elevationOne>=10*std::numbers::pi/180 && a.elevationTwo>=10*std::numbers::pi/180) || !std::isfinite(a.value)) throw std::runtime_error("Invalid imported VLBI observation");
            vlbi.push_back(a);++pairCounts[pair];
        }
    }
    std::sort(vlbi.begin(),vlbi.end(),[](const auto& a,const auto& b) {return a.elapsed<b.elapsed;});
    std::ofstream vlbiFile(o.output/"vlbi_observations.csv");vlbiFile<<std::setprecision(17)<<"elapsed_s,station_one,station_two,elevation_one_deg,elevation_two_deg,truth_delay_km,noise_m,observed_delay_km,sigma_m\n";
    for (const auto& a:vlbi) {const auto pair=pairs[a.pair];vlbiFile<<a.elapsed<<','<<names[pair[0]]<<','<<names[pair[1]]<<','<<a.elevationOne*180/std::numbers::pi<<','<<a.elevationTwo*180/std::numbers::pi<<','<<a.truth<<','<<a.noise*1000<<','<<a.value<<','<<o.vlbiSigma*1000<<'\n';}vlbiFile.close();
    std::vector<Eigen::Index> rangeRows(observations.size(),-1),dopplerRows(observations.size());
    Eigen::Index radioRows=0,rangeCount=0;
    std::vector<std::vector<Eigen::Index>> rangeGroups(o.passSchedule?static_cast<int>(o.days):0);
    for (std::size_t k=0;k<observations.size();++k) {if (observations[k].rangeActive) {rangeRows[k]=radioRows++;++rangeCount;if (o.passSchedule) rangeGroups[observations[k].pass].push_back(rangeRows[k]);}dopplerRows[k]=radioRows++;}
    rangeGroups.erase(std::remove_if(rangeGroups.begin(),rangeGroups.end(),[](const auto& g){return g.empty();}),rangeGroups.end());
    const Eigen::Index rows=radioRows+vlbi.size();Eigen::VectorXd observed(rows),variances(rows);
    for (std::size_t k=0;k<vlbi.size();++k) {observed[radioRows+k]=vlbi[k].value;variances[radioRows+k]=o.vlbiSigma*o.vlbiSigma;}
    for (std::size_t k=0;k<observations.size();++k) {const auto& a=observations[k];if (a.rangeActive) {observed[rangeRows[k]]=a.value[0];variances[rangeRows[k]]=o.rangeSigma*o.rangeSigma;}observed[dopplerRows[k]]=a.value[1];variances[dopplerRows[k]]=o.dopplerSigma*o.dopplerSigma;}
    const radio::PassRangeWhitening whitening(variances.cwiseSqrt(),rangeGroups,o.passSchedule?o.rangeBiasSigma:0);
    const Eigen::VectorXd fitVariances=o.passSchedule?Eigen::VectorXd::Ones(rows):variances;
    const auto weightedResidual=[&](const Eigen::VectorXd& r) {return whitening.vector(r);};
    std::ofstream obsFile(o.output/"observations.csv");obsFile<<std::setprecision(17)<<"elapsed_s,station,truth_range_km,truth_doppler_hz,range_noise_m,doppler_noise_hz,observed_range_km,observed_doppler_hz,pass_id,range_active,range_independent_noise_m,range_bias_m\n";
    for (const auto& a:observations) {obsFile<<a.elapsed<<','<<names[a.station]<<','<<a.truth[0]<<','<<a.truth[1]<<',';if (a.rangeActive) obsFile<<a.noise[0]*1000;obsFile<<','<<a.noise[1]<<',';if (a.rangeActive) obsFile<<a.value[0];obsFile<<','<<a.value[1]<<','<<a.pass<<','<<a.rangeActive<<','<<a.rangeIndependentNoise*1000<<','<<a.rangeBias*1000<<'\n';}obsFile.close();
    if (o.generateOnly) {std::cout<<"Generated fixed observation geometry and truth.\n";return 0;}
    std::cout<<observations.size()<<" tracking epochs ("<<rows<<" scalar observations).\n"<<std::flush;
    const Parameters offset=(Parameters()<<1,-.8,.5,1e-5,-8e-6,6e-6,0,0).finished(),seedState=trueInitial+o.initialGuessScale*offset;
    Parameters priorMean=trueInitial+offset;priorMean.tail<2>().setOnes();
    const Parameters priorSigma=(Parameters()<<2,2,2,1e-4,1e-4,1e-4,.2,.5).finished();
    std::vector<std::shared_ptr<Arc>> cache{nominal};
    const auto trajectory=[&](const Parameters& initial) {
        for (const auto& a:cache) if ((a->initial.array()==initial.array()).all()) return a;
        std::cout<<"  Nonlinear orbit propagation.\n"<<std::flush;if (initial[6]<=0 || initial[7]<=0) throw std::runtime_error("Negative force scale candidate");auto arc=std::make_shared<Arc>(t0,duration,initial,force,saturn,o.maxStep);
        cache.push_back(arc);if (cache.size()>3) cache.erase(cache.begin());return arc;
    };
    const auto physical=[&](const Eigen::VectorXd& z)->Parameters {return priorMean+priorSigma.cwiseProduct(z);};
    const auto evaluate=[&](const Eigen::VectorXd& z) {const auto arc=trajectory(physical(z));Eigen::VectorXd y(rows);for (std::size_t k=0;k<observations.size();++k) {const auto& a=observations[k];const auto predicted=value(links[a.station],t0+a.elapsed,o.count,*stations[a.station],*arc,ramp,referenceHz,clock);if (a.rangeActive) y[rangeRows[k]]=predicted[0];y[dopplerRows[k]]=predicted[1];}for (std::size_t k=0;k<vlbi.size();++k) {const auto& a=vlbi[k];const auto pair=pairs[a.pair];y[radioRows+k]=vlbiLinks[a.pair].solve(t0+a.elapsed,*arc,*stations[pair[0]],*stations[pair[1]]).differentialRangeKm;}return y;};
    // Affine sensitivity trajectories need no nonlinear force propagation at
    // these offsets. Large, symmetric radio stencils suppress GHz arithmetic
    // quantization, particularly in weak state/force parameter directions.
    Parameters differentiation;differentiation<<1000,1000,1000,1e-2,1e-2,1e-2,1000,1000;
    const auto raw=[&](const Eigen::VectorXd& z) {
        const auto arc=trajectory(physical(z));fd::filters::BatchLeastSquaresDriver::LinearizedProblem p;p.residuals=observed-evaluate(z);p.measurementVariances=fitVariances;p.designMatrix.resize(rows,8);
        const radio::TrajectoryFactory factory=[&](double,const Eigen::VectorXd& x) {return std::make_unique<AffineArc>(*arc,Parameters(x)-arc->initial);};
        for (std::size_t k=0;k<observations.size();++k) {const auto& a=observations[k];const radio::TwoWayMeasurement measurement(links[a.station],*stations[a.station],t0+a.elapsed,o.count,ramp,referenceHz,clock,factory,differentiation);const Eigen::MatrixXd H=measurement(t0,arc->initial).jacobian*priorSigma.asDiagonal();if (a.rangeActive) p.designMatrix.row(rangeRows[k])=H.row(0);p.designMatrix.row(dopplerRows[k])=H.row(1);}
        for (std::size_t k=0;k<vlbi.size();++k) {const auto& a=vlbi[k];const auto pair=pairs[a.pair];const radio::VlbiMeasurement m(vlbiLinks[a.pair],*stations[pair[0]],*stations[pair[1]],t0+a.elapsed,factory,differentiation);p.designMatrix.row(radioRows+k)=m(t0,arc->initial).jacobian*priorSigma.asDiagonal();}
        if (o.passSchedule) {p.residuals=whitening.vector(p.residuals);p.designMatrix=whitening.matrix(p.designMatrix);}
        return p;
    };
    const auto objective=[&](const Eigen::VectorXd& z,const Eigen::VectorXd& r) {return weightedResidual(r).squaredNorm()+z.squaredNorm();};
    const Eigen::VectorXd initialZ=(seedState-priorMean).cwiseQuotient(priorSigma);const Eigen::VectorXd prefit=observed-evaluate(initialZ);
    fd::filters::WLS wls;wls.setInitialState(initialZ,Matrix::Identity(),0);wls.setAprioriMean(Eigen::VectorXd::Zero(8));double damping=0;
    const auto linearize=[&](const Eigen::VectorXd& z) {
        const auto p=raw(z);const double old=p.residuals.cwiseQuotient(fitVariances.cwiseSqrt()).squaredNorm()+z.squaredNorm();
        auto undamped=wls;undamped.processBatchDiagonal(p.residuals,p.designMatrix,fitVariances);
        const Eigen::MatrixXd information=undamped.informationMatrix();
        for (int attempt=0;attempt<12;++attempt) {
            auto trial=p;if (damping>0) {trial.residuals.conservativeResize(rows+8);trial.residuals.tail<8>().setZero();trial.designMatrix.conservativeResize(rows+8,8);trial.designMatrix.bottomRows<8>()=std::sqrt(damping)*Matrix::Identity();trial.measurementVariances.conservativeResize(rows+8);trial.measurementVariances.tail<8>().setOnes();}
            auto probe=wls;probe.processBatchDiagonal(trial.residuals,trial.designMatrix,trial.measurementVariances);
            const Eigen::VectorXd delta=probe.lastStateCorrection();
            const double stepSigma=std::sqrt(std::max(0.,delta.dot(information*delta)));
            const double next=objective(probe.state(),observed-evaluate(probe.state()));
            // Counted Doppler at GHz frequencies has a finite precision floor
            // from differencing geometric delays. Only a step below 0.1
            // posterior sigma may use this 0.1% objective-resolution allowance.
            // Final convergence also requires a small UNdamped posterior step.
            const double allowance=stepSigma<.1 ? 1e-3*std::max(1.,old) : 1e-5;
            std::cout<<"  Trial: objective "<<old<<" -> "<<next<<", posterior step="<<stepSigma<<", damping="<<damping<<'\n'<<std::flush;
            if (next<=old+allowance) return trial;
            damping=damping==0?1:damping*10;
        }
        throw std::runtime_error("WLS failed to find decreasing objective");
    };
    fd::filters::BatchLeastSquaresDriver::ConvergenceCriteria criteria;criteria.maxIterations=o.iterations;criteria.minimumIterations=2;
    fd::filters::BatchLeastSquaresDriver driver(wls,linearize,criteria);
    driver.setConvergenceFunction([&](const auto& it,const auto& c) {
        const auto p=raw(it.stateAfter);auto check=wls;
        check.processBatchDiagonal(p.residuals,p.designMatrix,fitVariances);
        const Eigen::VectorXd delta=check.lastStateCorrection();
        const double posteriorStep=std::sqrt(std::max(0.,delta.dot(check.informationMatrix()*delta)));
        std::cout<<"  Undamped posterior stationarity="<<posteriorStep<<" sigma\n"<<std::flush;
        damping*=.1;
        return it.iteration>=c.minimumIterations && posteriorStep<.1;
    });
    std::ofstream iterations(o.output/"iterations.csv");iterations<<std::setprecision(17)<<"iteration,weighted_prefit_rms,weighted_postfit_rms,normalized_correction_norm,regularized_objective,damping\n";
    const auto result=driver.solve([&](const auto& it) {const Eigen::VectorXd r=observed-evaluate(it.stateAfter);const double rms=weightedResidual(r).norm()/std::sqrt(rows);iterations<<it.iteration<<','<<it.weightedRms<<','<<rms<<','<<it.correctionNorm<<','<<objective(it.stateAfter,r)<<','<<damping<<'\n';iterations.flush();std::cout<<"Iteration "<<it.iteration<<": weighted RMS="<<rms<<", correction="<<it.correctionNorm<<'\n'<<std::flush;});iterations.close();
    const auto finalProblem=raw(wls.state());auto finalSolver=wls;finalSolver.processBatchDiagonal(finalProblem.residuals,finalProblem.designMatrix,fitVariances);
    const Eigen::VectorXd finalCorrection=finalSolver.lastStateCorrection();
    const double undampedCorrection=finalCorrection.norm();
    const double undampedPosteriorStep=std::sqrt(std::max(0.,finalCorrection.dot(finalSolver.informationMatrix()*finalCorrection)));
    const bool converged=result.converged() && undampedPosteriorStep<.1;
    const Matrix covariance=priorSigma.asDiagonal()*finalSolver.covariance()*priorSigma.asDiagonal();const Parameters estimate=physical(wls.state());const auto fitted=trajectory(estimate),seedArc=trajectory(seedState);
    const Eigen::VectorXd postfit=observed-evaluate(wls.state());
    // Verify finite-difference radio sensitivity against a second stencil.
    double stencilRelative=0,stencilPosteriorNoise=0;
    const Eigen::MatrixXd posteriorRoot=finalSolver.covariance().llt().matrixL();
    for (std::size_t k:{std::size_t(0),observations.size()/2,observations.size()-1}) {
        const auto& a=observations[k];const radio::TrajectoryFactory factory=[&](double,const Eigen::VectorXd& x) {return std::make_unique<AffineArc>(*fitted,Parameters(x)-estimate);};
        const radio::TwoWayMeasurement full(links[a.station],*stations[a.station],t0+a.elapsed,o.count,ramp,referenceHz,clock,factory,differentiation),half(links[a.station],*stations[a.station],t0+a.elapsed,o.count,ramp,referenceHz,clock,factory,differentiation/2);
        const Eigen::MatrixXd A=full(t0,estimate).jacobian*priorSigma.asDiagonal(),B=half(t0,estimate).jacobian*priorSigma.asDiagonal();
        for (int row=0;row<2;++row) {
            if (row==0 && !a.rangeActive) continue;
            stencilRelative=std::max(stencilRelative,(A.row(row)-B.row(row)).norm()/std::max(1e-20,A.row(row).norm()));
            stencilPosteriorNoise=std::max(stencilPosteriorNoise,((A.row(row)-B.row(row))*posteriorRoot).norm()/std::sqrt(variances[row==0?rangeRows[k]:dopplerRows[k]]));
        }
    }
    double vlbiStencilNoise=0;
    for (std::size_t k=0;k<vlbi.size();k+=std::max(std::size_t(1),vlbi.size()/3)) {
        const auto& a=vlbi[k];const auto pair=pairs[a.pair];const radio::TrajectoryFactory factory=[&](double,const Eigen::VectorXd& x) {return std::make_unique<AffineArc>(*fitted,Parameters(x)-estimate);};
        const radio::VlbiMeasurement full(vlbiLinks[a.pair],*stations[pair[0]],*stations[pair[1]],t0+a.elapsed,factory,differentiation),half(vlbiLinks[a.pair],*stations[pair[0]],*stations[pair[1]],t0+a.elapsed,factory,differentiation/2);
        vlbiStencilNoise=std::max(vlbiStencilNoise,((full(t0,estimate).jacobian-half(t0,estimate).jacobian)*priorSigma.asDiagonal()*posteriorRoot).norm()/o.vlbiSigma);
    }
    if (vlbiStencilNoise>.01) throw std::runtime_error("VLBI stencil discrepancy exceeds 1% of noise over posterior uncertainty");
    if (stencilRelative>.03) throw std::runtime_error("Radio derivative stencil check exceeds 3%");
    if (stencilPosteriorNoise>.01) throw std::runtime_error("Radio stencil discrepancy exceeds 1% of measurement noise over posterior uncertainty");
    std::cout<<"Validating final orbit with half the maximum integration step.\n"<<std::flush;
    const Arc refined(t0,duration,estimate,force,saturn,o.maxStep/2);double numericalMaximum=0,rangeNumerical=0,dopplerNumerical=0;
    std::ofstream residualFile(o.output/"residuals.csv");residualFile<<std::setprecision(17)<<"elapsed_s,station,prefit_range_m,postfit_range_m,prefit_doppler_hz,postfit_doppler_hz\n";
    for (std::size_t k=0;k<observations.size();++k) {const auto& a=observations[k];residualFile<<a.elapsed<<','<<names[a.station]<<',';if (a.rangeActive) residualFile<<prefit[rangeRows[k]]*1000;residualFile<<',';if (a.rangeActive) residualFile<<postfit[rangeRows[k]]*1000;residualFile<<','<<prefit[dopplerRows[k]]<<','<<postfit[dopplerRows[k]]<<'\n';if (k%50==0 || k==observations.size()-1) {const auto refinedValue=value(links[a.station],t0+a.elapsed,o.count,*stations[a.station],refined,ramp,referenceHz,clock);const auto actual=value(links[a.station],t0+a.elapsed,o.count,*stations[a.station],*fitted,ramp,referenceHz,clock);if (a.rangeActive) rangeNumerical=std::max(rangeNumerical,std::abs(refinedValue[0]-actual[0])*1000);dopplerNumerical=std::max(dopplerNumerical,std::abs(refinedValue[1]-actual[1]));}}residualFile.close();
    double vlbiNumerical=0,vlbiSum=0;
    std::ofstream vlbiResidual(o.output/"vlbi_residuals.csv");vlbiResidual<<std::setprecision(17)<<"elapsed_s,station_one,station_two,prefit_m,postfit_m\n";
    for (std::size_t k=0;k<vlbi.size();++k) {const auto& a=vlbi[k];const auto pair=pairs[a.pair];vlbiResidual<<a.elapsed<<','<<names[pair[0]]<<','<<names[pair[1]]<<','<<prefit[radioRows+k]*1000<<','<<postfit[radioRows+k]*1000<<'\n';vlbiSum+=postfit[radioRows+k]*postfit[radioRows+k];
        const double refinedValue=vlbiLinks[a.pair].solve(t0+a.elapsed,refined,*stations[pair[0]],*stations[pair[1]]).differentialRangeKm,actual=vlbiLinks[a.pair].solve(t0+a.elapsed,*fitted,*stations[pair[0]],*stations[pair[1]]).differentialRangeKm;
        vlbiNumerical=std::max(vlbiNumerical,std::abs(refinedValue-actual)*1000);
    }vlbiResidual.close();
    if (vlbiNumerical>o.vlbiSigma*1000*.01) throw std::runtime_error("VLBI refined propagation exceeds 1% of noise");
    std::ofstream initialFile(o.output/"initial_estimate.csv");initialFile<<std::setprecision(17)<<"component,unit,truth,initial_guess,estimate,posterior_sigma\n";
    const std::array<const char*,8> componentNames{"x","y","z","vx","vy","vz","C_s","lambda"};
    for (int j=0;j<8;++j) initialFile<<componentNames[j]<<','<<(j<3?"km":j<6?"km/s":"dimensionless")<<','<<trueInitial[j]<<','<<seedState[j]<<','<<estimate[j]<<','<<std::sqrt(covariance(j,j))<<'\n';initialFile.close();
    std::ofstream orbit(o.output/"trajectory.csv");orbit<<std::setprecision(17)<<"elapsed_s,truth_x_km,truth_y_km,truth_z_km,nominal_x_km,nominal_y_km,nominal_z_km,fitted_x_km,fitted_y_km,fitted_z_km,nominal_position_error_m,seed_position_error_m,fitted_position_error_m,fitted_velocity_error_mm_s,radial_error_m,along_error_m,cross_error_m,formal_position_sigma_m,fitted_vx_km_s,fitted_vy_km_s,fitted_vz_km_s\n";
    double nominalSum=0,seedSum=0,fittedSum=0,endpoint=0,nominalEndpoint=0;std::size_t samples=0;
    std::ofstream epochCov(o.output/"epoch_error_covariance.csv");epochCov<<std::setprecision(17)<<"elapsed_s,dx_m,dy_m,dz_m,dvx_mm_s,dvy_mm_s,dvz_mm_s";
    for (int i=0;i<6;++i) for (int j=0;j<6;++j) epochCov<<",cov_"<<i<<j;epochCov<<'\n';
    for (std::size_t k=0;;++k) {const double elapsed=std::min(k*o.cadence,duration),t=t0+elapsed;const State target=truth.relative(t),x=fitted->relative(t),base=nominal->relative(t),seed=seedArc->relative(t);const State error=x-target;const double ne=(base.head<3>()-target.head<3>()).norm()*1000,se=(seed.head<3>()-target.head<3>()).norm()*1000,fe=error.head<3>().norm()*1000;nominalSum+=ne*ne;seedSum+=se*se;fittedSum+=fe*fe;endpoint=fe;nominalEndpoint=ne;++samples;
        const Eigen::Vector3d radial=target.head<3>().normalized(),normal=target.head<3>().cross(target.tail<3>()).normalized(),along=normal.cross(radial);const Sensitivity phi=fitted->transition(t);const Eigen::Matrix<double,6,6> P=phi*covariance*phi.transpose();numericalMaximum=std::max(numericalMaximum,(x.head<3>()-refined.relative(t).head<3>()).norm()*1000);
        if (std::fmod(elapsed,86400)==0 || elapsed==duration) {
            epochCov<<elapsed;for (int i=0;i<6;++i) epochCov<<','<<error[i]*(i<3?1000:1e6);
            for (int i=0;i<6;++i) for (int j=0;j<6;++j) epochCov<<','<<P(i,j)*(i<3?1000:1e6)*(j<3?1000:1e6);epochCov<<'\n';
        }
        orbit<<elapsed;for (const State& state:{target,base,x}) for (int j=0;j<3;++j) orbit<<','<<state[j];orbit<<','<<ne<<','<<se<<','<<fe<<','<<error.tail<3>().norm()*1e6<<','<<error.head<3>().dot(radial)*1000<<','<<error.head<3>().dot(along)*1000<<','<<error.head<3>().dot(normal)*1000<<','<<std::sqrt(P.topLeftCorner<3,3>().trace())*1000;for (int j=3;j<6;++j) orbit<<','<<x[j];orbit<<'\n';if (elapsed==duration) break;
    }orbit.close();
    if (numericalMaximum>.1 || rangeNumerical>.1 || dopplerNumerical>o.dopplerSigma*.1) throw std::runtime_error("Numerical convergence exceeds baseline thresholds");
    std::ofstream cov(o.output/"initial_covariance.csv");cov<<std::setprecision(17);for (int i=0;i<8;++i) {for (int j=0;j<8;++j) cov<<(j?",":"")<<covariance(i,j);cov<<'\n';}cov.close();
    Eigen::MatrixXd weighted=finalProblem.designMatrix;for (Eigen::Index i=0;i<rows;++i) weighted.row(i)/=std::sqrt(fitVariances[i]);Eigen::JacobiSVD<Eigen::MatrixXd> svd(weighted);const auto singular=svd.singularValues();
    const auto rmsComponent=[&](const Eigen::VectorXd& r,int j) {double sum=0;for (std::size_t k=0;k<observations.size();++k) {if (j==0 && !observations[k].rangeActive) continue;const auto row=j==0?rangeRows[k]:dopplerRows[k];sum+=r[row]*r[row];}return std::sqrt(sum/(j==0?rangeCount:observations.size()));};
    const double runtime=std::chrono::duration<double>(std::chrono::steady_clock::now()-wall).count();
    std::ofstream summary(o.output/"summary.json");summary<<std::setprecision(17)<<"{\n\"truth\":\""<<o.truth<<"\",\n\"start_utc\":\""<<o.utc<<"\",\n\"duration_days\":"<<o.days<<",\n\"cadence_s\":"<<o.cadence<<",\n\"count_s\":"<<o.count<<",\n\"mass_kg\":"<<o.mass<<",\n\"area_m2\":"<<o.area<<",\n\"cr\":"<<o.cr<<",\n\"srp_prior_sigma\":0.2,\n\"thermal_prior_sigma\":0.5,\n\"truth_source\":\"results/cassini_two_way_wls/pass_schedule/truth/trajectory.csv:multibody_j2_srp_thermal\",\n\"gravity_mask\":2,\n\"thermal_reference_vector_km_s2\":";jsonVector(summary,o.recoil);
    summary<<",\n\"pass_schedule\":"<<(o.passSchedule?"true":"false")<<",\n\"pass_count\":"<<(o.passSchedule?static_cast<int>(o.days):0)<<",\n\"pass_target_hours\":"<<o.passHours<<",\n\"range_observations\":"<<rangeCount<<",\n\"doppler_observations\":"<<observations.size()<<",\n\"range_cadence_s\":"<<o.rangeCadence<<",\n\"range_independent_sigma_m\":"<<o.rangeSigma*1000<<",\n\"range_pass_bias_sigma_m\":"<<(o.passSchedule?o.rangeBiasSigma*1000:0);
    summary<<",\n\"initial_guess_scale\":"<<o.initialGuessScale;
    summary<<",\n\"vlbi_observations\":"<<vlbi.size()<<",\n\"vlbi_sigma_m\":"<<o.vlbiSigma*1000<<",\n\"vlbi_noise_seed\":"<<(o.observationInput.empty()?std::to_string(static_cast<std::uint64_t>(o.seed)+1000003):"null")<<",\n\"vlbi_pair_counts\":["<<pairCounts[0]<<','<<pairCounts[1]<<','<<pairCounts[2]<<"],\n\"postfit_vlbi_rms_m\":"<<(vlbi.empty()?0:std::sqrt(vlbiSum/vlbi.size())*1000)<<",\n\"vlbi_stencil_posterior_noise_fraction\":"<<vlbiStencilNoise<<",\n\"refined_vlbi_maximum_m\":"<<vlbiNumerical;
    summary<<",\n\"noise_seed\":"<<(o.observationInput.empty()?std::to_string(o.seed):"null")<<",\n\"noise_enabled\":"<<(o.noise?"true":"false")<<",\n\"range_sigma_m\":"<<std::hypot(o.rangeSigma,o.passSchedule?o.rangeBiasSigma:0)*1000<<",\n\"doppler_sigma_hz\":"<<o.dopplerSigma<<",\n\"candidate_epochs\":"<<candidates<<",\n\"tracking_epochs\":"<<observations.size()<<",\n\"scalar_observations\":"<<rows<<",\n\"station_counts\":["<<stationCounts[0]<<','<<stationCounts[1]<<','<<stationCounts[2]<<"],\n\"estimated_parameters\":8,\n\"iterations\":"<<result.iterationCount()<<",\n\"converged\":"<<(converged?"true":"false")<<",\n\"undamped_correction_norm\":"<<undampedCorrection<<",\n\"undamped_posterior_step_sigma\":"<<undampedPosteriorStep<<",\n\"small_step_objective_relative_allowance\":0.001,\n\"prefit_range_rms_m\":"<<rmsComponent(prefit,0)*1000<<",\n\"postfit_range_rms_m\":"<<rmsComponent(postfit,0)*1000<<",\n\"prefit_doppler_rms_hz\":"<<rmsComponent(prefit,1)<<",\n\"postfit_doppler_rms_hz\":"<<rmsComponent(postfit,1)<<",\n\"weighted_postfit_rms\":"<<weightedResidual(postfit).norm()/std::sqrt(rows)<<",\n\"reduced_chi_square\":"<<weightedResidual(postfit).squaredNorm()/(rows-8)<<",\n\"nominal_position_rms_m\":"<<std::sqrt(nominalSum/samples)<<",\n\"seed_position_rms_m\":"<<std::sqrt(seedSum/samples)<<",\n\"fitted_position_rms_m\":"<<std::sqrt(fittedSum/samples)<<",\n\"nominal_endpoint_error_m\":"<<nominalEndpoint<<",\n\"fitted_endpoint_error_m\":"<<endpoint<<",\n\"nominal_refined_maximum_m\":"<<numericalMaximum<<",\n\"nominal_refined_range_maximum_m\":"<<rangeNumerical<<",\n\"nominal_refined_doppler_maximum_hz\":"<<dopplerNumerical<<",\n\"radio_stencil_relative_difference\":"<<stencilRelative<<",\n\"measurement_design_condition_number\":"<<singular[0]/singular[7]<<",\n\"measurement_design_singular_values\":";jsonVector(summary,singular);
    summary<<",\n\"radio_stencil_posterior_noise_fraction\":"<<stencilPosteriorNoise<<",\n\"radio_differentiation_steps\":";jsonVector(summary,differentiation);summary<<",\n\"initial_prior_mean\":";jsonVector(summary,priorMean);summary<<",\n\"initial_truth_state\":";jsonVector(summary,trueInitial);summary<<",\n\"initial_seed_state\":";jsonVector(summary,seedState);summary<<",\n\"initial_estimated_state\":";jsonVector(summary,estimate);summary<<",\n\"initial_prior_sigma\":";jsonVector(summary,priorSigma);summary<<",\n\"initial_posterior_sigma\":";jsonVector(summary,covariance.diagonal().cwiseSqrt());summary<<",\n\"runtime_s\":"<<runtime<<"\n}\n";summary.close();
    std::cout<<"Finished: converged="<<converged<<", 3D RMS="<<std::sqrt(fittedSum/samples)<<" m; numerical maximum="<<numericalMaximum<<" m.\n";
    return converged?0:2;
 } catch (const std::exception& e) {std::cerr<<e.what()<<'\n';return 1;}
}
