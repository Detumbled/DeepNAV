#include "../cassini_trajectory/SavedSimulation.hpp"
#include "opnav/synth/SyntheticImageGenerator.hpp"
#include "opnav/estimation/BodyEstimator.hpp"
#include "opnav/estimation/SpiceBodyModel.hpp"
#include "opnav/core/ApparentDirection.hpp"
#include "dynamics/SpiceKernelSet.hpp"
#include "dynamics/SpiceStateProvider.hpp"
#include "simulation/CassiniRadioKernels.hpp"
#include <Eigen/Geometry>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <numbers>

namespace fs=std::filesystem;
using namespace fd::opnav;
std::ofstream output(const fs::path& path){std::ofstream out(path);out.exceptions(std::ios::badbit|std::ios::failbit);out<<std::setprecision(17);return out;}
struct Scenario {
    std::string name;bool psf=false,shape=false,surface=false,calibration=false,attitude=false,exposure=false,ellipsoidFit=false,integrate=false;
};
int main(int argc,char** argv){try{
    const fs::path root=DEEPNAV_SOURCE_DIR;fs::path directory=root/"results/cassini_synthetic_optical/realism";int trials=24;
    for(int i=1;i<argc;++i){const std::string key=argv[i];if(key=="--help"){std::cout<<"cassini_optical_realism [--output DIR] [--trials N]\n";return 0;}
        if(++i==argc)throw std::invalid_argument("Missing argument.");if(key=="--output")directory=fs::absolute(argv[i]);
        else if(key=="--trials")trials=std::stoi(argv[i]);else throw std::invalid_argument("Unknown argument.");}
    if(trials<1 || trials>1000)throw std::invalid_argument("Invalid Monte Carlo trial count.");
    od::SpiceErrorModeGuard errors("RETURN","NONE");auto paths=fd::simulation::cassiniRadioKernels(root/"Kernels");
    for(const char* file:{"Cassini/fk/cas_v40.tf","Cassini/sclk/cas00172.tsc","Cassini/ck/04256_04292py_as_flown.bc",
        "Cassini/ck/04292_04320py_as_flown.bc","Cassini/ck/04281_04286ra.bc","Cassini/ck/04286_04291ra.bc",
        "Cassini/ck/04291_04296ra.bc","Cassini/ck/04296_04301ra.bc","Cassini/ck/04293_04298ra.bc"})paths.push_back(root/"Kernels"/file);
    const fd::dynamics::SpiceKernelSet kernels(paths);const CassiniAttitude forceAttitude(paths);
    const Sat427Gravity gravity(root/"Kernels/Cassini/spk/200128RU_SCPSE_04251_04303.bsp.lbl");
    double recoilEpoch;str2et_c("2005-08-13T00:00:00",&recoilEpoch);double radii[3];SpiceInt count;bodvrd_c("SATURN","RADII",3,&count,radii);
    od::throwIfSpiceFailed("Realism force setup");
    const fd::simulation::CassiniForceModel force(gravity,forceAttitude,fd::simulation::cassiniPerturbers(gravity),2,3220,12,1.42,
        {-2.725e-13,-2.725e-13,-5.232237415805135e-12},recoilEpoch,radii[0]);
    const fd::dynamics::SpiceStateProvider saturn("SATURN"),target("ENCELADUS");
    const fs::path truthPath=root/"results/cassini_two_way_wls/pass_schedule/truth/trajectory.csv";
    const fd::simulation::SavedSimulation truth(truthPath,force,saturn);
    double epoch;str2et_c("2004-10-10T00:00:00",&epoch);od::throwIfSpiceFailed("Realism epoch");
    const auto direction=core::computeApparentDirection(TdbEpoch(epoch),truth.stateAt(TdbEpoch(epoch)),target);
    const Eigen::Vector3d z=direction.unitDirection,x=z.unitOrthogonal(),y=z.cross(x);Eigen::Matrix3d base;
    base.row(0)=x.transpose();base.row(1)=y.transpose();base.row(2)=z.transpose();const double scale=2003.44/.012;
    // Off-axis target makes camera-calibration errors observable in this study.
    base=Eigen::AngleAxisd(5/scale,Eigen::Vector3d::UnitX()).toRotationMatrix()*
         Eigen::AngleAxisd(-10/scale,Eigen::Vector3d::UnitY()).toRotationMatrix()*base;
    const CameraIntrinsics intrinsics{2003.44,{31.23,30.74},-Eigen::Matrix2d::Identity()/.012};
    const CameraModel nominalCamera(intrinsics);
    const std::vector<Scenario> scenarios={
        {"matched"},{"psf",true},{"shape_sphere",false,true},
        {"shape_ellipsoid",false,true,false,false,false,false,true},
        {"surface",false,false,true},{"calibration",false,false,false,true},
        {"attitude",false,false,false,false,true},
        {"exposure_midpoint",false,false,false,false,false,true},
        {"exposure_integrated",false,false,false,false,false,true,false,true},
        {"combined_sphere",true,true,true,true,true,true,false,true},
        {"combined_ellipsoid",true,true,true,true,true,true,true,true}};
    fs::create_directories(directory);auto results=output(directory/"cases.csv"),mc=output(directory/"trials.csv");
    const std::string header="case,trial,status,fit_checks_pass,model_mismatch,center_sample,center_line,extract_dx,extract_dy,nav_dx,nav_dy,chi2,formal_ss,formal_sl,formal_ll,effective_ss,effective_sl,effective_ll,gaussian_error\n";
    results<<header;mc<<header;
    auto config=output(directory/"config.json");config<<"{\n\"utc\":\"2004-10-10T00:00:00\",\n\"truth_csv\":"<<std::quoted(truthPath.string())
        <<",\n\"trials\":"<<trials<<",\n\"pixel_scale_per_radian\":"<<scale
        <<",\n\"nominal_psf_sigma\":0.8,\n\"truth_psf_sample_line\":[0.88,0.72],\n\"albedo_contrast\":0.12,\n\"lommel_seeliger_mix\":0.15,"
        <<"\n\"calibration\":{\"relative_focal_error\":0.0004,\"principal_error\":[0.12,-0.08],\"distortion_k1\":50000},"
        <<"\n\"attitude\":{\"equivalent_static_pixel_error\":[0.25,0.15],\"jitter_amplitude_pixels\":0.08,\"nominal_sigma_pixels\":0.25},"
        <<"\n\"nominal_detector\":{\"electrons_per_dn\":10,\"read_noise_electrons\":3,\"bias_dn\":100},"
        <<"\n\"long_exposure_seconds\":4,\n\"nuisance_draws\":\"fixed declared offsets; Monte Carlo varies detector noise only\","
        <<"\n\"calibration_parameters\":\"illustrative synthetic stress settings, not flight-calibrated Cassini values\",\n\"loaded_kernels\":[";
    for(std::size_t i=0;i<paths.size();++i)config<<(i?",":"")<<std::quoted(paths[i].string());config<<"]\n}\n";
    for(const auto& scenario:scenarios){
        synth::SyntheticImageOptions options;options.size={64,64};options.subpixels=8;options.psfSigmaPixels=scenario.psf?.88:.8;
        options.psfSigmaLinePixels=scenario.psf?.72:.8;options.usePckEllipsoid=scenario.shape;
        options.albedoContrast=scenario.surface?.12:0;options.lommelSeeligerMix=scenario.surface?.15:0;
        options.exposureSeconds=scenario.exposure?4:.06;options.exposureSamples=scenario.exposure?21:5;
        options.referenceRangeKm=direction.lightTime.geometricLineOfSightKm.norm();options.referenceExposureSeconds=options.exposureSeconds;
        CameraIntrinsics actual=intrinsics;if(scenario.calibration){actual.focalLength*=1.0004;actual.principalPoint.sample+=.12;actual.principalPoint.line-=.08;}
        OpenCVDistortionCoefficients distortion;distortion.K1=scenario.calibration?50000:0;
        const CameraModel actualCamera(actual,std::make_shared<OpenCVDistortionModel>(distortion));
        const auto nominalRotation=[=](double t)->Eigen::Matrix3d {
            const double dt=t-epoch;
            return Eigen::AngleAxisd((scenario.exposure?.20:0)*dt/scale,Eigen::Vector3d::UnitX()).toRotationMatrix()*
                   Eigen::AngleAxisd(-(scenario.exposure?.45:0)*dt/scale,Eigen::Vector3d::UnitY()).toRotationMatrix()*base;
        };
        const auto actualAttitude=[=](TdbEpoch t){
            const double dt=t.secondsPastJ2000-epoch;
            const double jitter=scenario.attitude?.08*std::sin(2*std::numbers::pi*dt/.15):0;
            const Eigen::Matrix3d offset=Eigen::AngleAxisd((scenario.attitude?.15:0)/scale,Eigen::Vector3d::UnitX()).toRotationMatrix()*
                Eigen::AngleAxisd(-((scenario.attitude?.25:0)+jitter)/scale,Eigen::Vector3d::UnitY()).toRotationMatrix();
            return CameraAttitude(offset*nominalRotation(t.secondsPastJ2000));
        };
        synth::OpticalDetectorNoise noise;noise.electronsPerDn=10;noise.backgroundElectrons=50;noise.biasDn=100;noise.seed=7100;
        const synth::SyntheticImageGenerator generator(actualCamera,actualAttitude,options);
        const auto frame=generator.createImage("ENCELADUS",TdbEpoch(epoch),noise,truth,"SIMULATED_CASSINI");
        OpticalImage observation;observation.body="ENCELADUS";observation.time=TdbEpoch(epoch);observation.exposureSeconds=options.exposureSeconds;
        observation.attitude=CameraAttitude(nominalRotation(epoch));observation.dn=frame.dn;observation.mask=frame.mask;
        if(scenario.attitude)observation.attitudeCovariance=Eigen::Matrix3d::Identity()*std::pow(.25/scale,2);
        const auto nominalPixel=[&](double t){const auto d=core::computeApparentDirection(TdbEpoch(t),truth.stateAt(TdbEpoch(t)),target);
            return nominalCamera.project(d.vectorKm,CameraAttitude(nominalRotation(t)));};
        const auto before=nominalPixel(epoch-options.exposureSeconds/2),after=nominalPixel(epoch+options.exposureSeconds/2),predicted=nominalPixel(epoch);
        observation.smearPixels={after.sample-before.sample,after.line-before.line};
        CameraCalibration calibration;calibration.intrinsics=intrinsics;
        if(scenario.calibration){calibration.principalPointCovariance=Eigen::Matrix2d::Identity()*.01;calibration.relativeFocalSigma=.0005;}
        BodyModel body;if(scenario.ellipsoidFit){body=makeSpiceBodyModel("ENCELADUS",direction.lightTime.emissionEpoch,
            observation.attitude,intrinsics);observation.bodyToView=body.bodyToView;}
        BodyFitOptions fitOptions;fitOptions.exposureSamples=scenario.integrate?9:1;
        const BodyEstimator estimator(BodyRenderer(body,16),calibration,fitOptions);
        const auto write=[&](std::ofstream& file,int trial,const BodyFitResult& fit){
            const double nan=std::numeric_limits<double>::quiet_NaN();const Eigen::Matrix2d effective=fit.measurement?fit.measurement->covariance:Eigen::Matrix2d::Constant(nan);
            const image::Image inferred=(observation.dn.array()-calibration.biasDn).max(0)/calibration.electronsPerDn+
                std::pow(calibration.readNoiseElectrons/calibration.electronsPerDn,2)+1./12;
            const auto gaussian=image::fitBrightestCircularGaussian(observation.dn,inferred,12,{},&observation.mask);
            const double gaussianError=gaussian.measurement?std::hypot(gaussian.model.center.sample-frame.truth.geometricCenter.sample,gaussian.model.center.line-frame.truth.geometricCenter.line):nan;
            file<<scenario.name<<','<<trial<<','<<int(fit.status)<<','<<fit.passesFitChecks()<<','<<fit.modelMismatch<<','<<fit.parameters.center.sample<<','<<fit.parameters.center.line<<','
                <<fit.parameters.center.sample-frame.truth.geometricCenter.sample<<','<<fit.parameters.center.line-frame.truth.geometricCenter.line<<','
                <<fit.parameters.center.sample-predicted.sample<<','<<fit.parameters.center.line-predicted.line<<','<<fit.reducedChiSquared<<','
                <<fit.formalCenterCovariance(0,0)<<','<<fit.formalCenterCovariance(0,1)<<','<<fit.formalCenterCovariance(1,1)<<','
                <<effective(0,0)<<','<<effective(0,1)<<','<<effective(1,1)<<','<<gaussianError<<'\n';
        };
        const auto fit=estimator.fit(observation);write(results,0,fit);results.flush();
        auto pixels=output(directory/(scenario.name+".csv"));pixels<<"sample,line,dn,source_electrons,variance_estimated,model_dn,residual_dn,truth_sample,truth_line,nominal_sample,nominal_line\n";
        for(int l=0;l<fit.modelDn.rows();++l)for(int s=0;s<fit.modelDn.cols();++s){const int x=int(fit.origin.sample)+s,y=int(fit.origin.line)+l;
            pixels<<x<<','<<y<<','<<observation.dn(y,x)<<','<<frame.expectedSourceElectrons(y,x)<<','<<fit.varianceDn2(l,s)<<','<<fit.modelDn(l,s)<<','<<fit.residualDn(l,s)<<','
                  <<frame.truth.geometricCenter.sample<<','<<frame.truth.geometricCenter.line<<','<<predicted.sample<<','<<predicted.line<<'\n';}
        std::cout<<scenario.name<<": center error "<<std::hypot(fit.parameters.center.sample-frame.truth.geometricCenter.sample,fit.parameters.center.line-frame.truth.geometricCenter.line)
                 <<" px, navigation error "<<std::hypot(fit.parameters.center.sample-predicted.sample,fit.parameters.center.line-predicted.line)<<" px, chi2="<<fit.reducedChiSquared<<", fit checks="<<fit.passesFitChecks()<<'\n'<<std::flush;
        if(scenario.name=="matched" || scenario.name=="combined_ellipsoid")for(int i=0;i<trials;++i){noise.seed=9000+i;
            const auto sampled=synth::sampleDetector(frame.expectedSourceElectrons,noise);observation.dn=sampled.dn;observation.mask=sampled.mask;
            const auto result=estimator.fit(observation);write(mc,i,result);mc.flush();}
    }
    std::cout<<"Saved C++ fits and robustness cases in "<<directory<<'\n';
}catch(const std::exception& e){std::cerr<<"ERROR: "<<e.what()<<'\n';return 1;}}
