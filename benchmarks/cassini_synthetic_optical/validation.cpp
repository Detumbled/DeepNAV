#include "../cassini_trajectory/SavedSimulation.hpp"
#include "opnav/synth/SyntheticImageGenerator.hpp"
#include "opnav/estimation/BodyEstimator.hpp"
#include "opnav/estimation/EstimatorRejection.hpp"
#include "opnav/estimation/SpiceBodyModel.hpp"
#include "opnav/core/ApparentDirection.hpp"
#include "dynamics/SpiceKernelSet.hpp"
#include "dynamics/SpiceStateProvider.hpp"
#include "simulation/CassiniRadioKernels.hpp"
#include <Eigen/Geometry>
#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numbers>
#include <random>

namespace fs = std::filesystem;
using namespace fd::opnav;
constexpr double missingValue = std::numeric_limits<double>::quiet_NaN();

std::ofstream output(const fs::path& path) {
    std::ofstream stream(path);
    stream.exceptions(std::ios::badbit | std::ios::failbit);
    stream << std::setprecision(17);
    return stream;
}

struct Nuisance {
    Eigen::Vector2d principal{Eigen::Vector2d::Zero()}, alignment{Eigen::Vector2d::Zero()};
    double focal=0, k1=0, psfX=.8, psfY=.8, albedo=0, reflectance=0;
    double jitter=0, jitterPhase=0, gain=10, read=3, bias=100;
};

// Separate, explicit distributions make the covariance experiment reproducible.
Nuisance drawNuisance(std::mt19937_64& random, double alignmentSigma) {
    std::normal_distribution<double> normal;
    std::uniform_real_distribution<double> uniform;
    Nuisance n;
    n.principal={.1*normal(random),.1*normal(random)};
    n.alignment={alignmentSigma*normal(random),alignmentSigma*normal(random)};
    n.focal=.0005*normal(random); n.k1=30000*normal(random);
    const double width=.04*normal(random), ellipticity=.08*normal(random);
    n.psfX=.8*std::exp(width+ellipticity); n.psfY=.8*std::exp(width-ellipticity);
    n.albedo=.03+.17*uniform(random); n.reflectance=.05+.20*uniform(random);
    n.jitter=std::abs(.08*normal(random)); n.jitterPhase=2*std::numbers::pi*uniform(random);
    n.gain=10*std::exp(.01*normal(random)); n.read=3*std::exp(.10*normal(random));
    n.bias=100+.1*normal(random);
    return n;
}

struct Scene {
    int id, day, draw, track;
    std::string experiment;
    std::uint64_t detectorSeed;
    double epoch, phaseDegrees, radiusPixels, budget;
    Eigen::Matrix3d nominalRotation;
    CameraIntrinsics intrinsics;
    Eigen::Vector2d pointing;
    Nuisance nuisance;
    std::optional<Eigen::Vector3d> controlledSun;
    core::ApparentDirection direction{{TdbEpoch(0),TdbEpoch(0)}};
};

struct Observation {
    synth::SyntheticOpticalFrame frame;
    OpticalImage image;
    CameraCalibration calibration;
    PixelCoordinates nominalPixel;
};

Observation render(const Scene& scene, bool perturbed, int subpixels,
                   const fd::dynamics::StateProvider& truth) {
    const Nuisance n=perturbed ? scene.nuisance : Nuisance{};
    CameraIntrinsics actual=scene.intrinsics;
    actual.focalLength*=1+n.focal;
    actual.principalPoint.sample+=n.principal.x(); actual.principalPoint.line+=n.principal.y();
    OpenCVDistortionCoefficients distortion; distortion.K1=n.k1;
    const CameraModel camera(actual,std::make_shared<OpenCVDistortionModel>(distortion));
    const double scale=scene.intrinsics.focalLength/.012;
    const Eigen::Vector2d pointing=perturbed ? scene.pointing : Eigen::Vector2d::Zero();
    const auto offset=[=](double dt)->Eigen::Matrix3d {
        const double jitter=n.jitter*std::sin(2*std::numbers::pi*dt/.15+n.jitterPhase);
        return Eigen::AngleAxisd(pointing.y()/scale,Eigen::Vector3d::UnitX()).toRotationMatrix()*
            Eigen::AngleAxisd(-(pointing.x()+jitter)/scale,Eigen::Vector3d::UnitY()).toRotationMatrix();
    };
    const auto attitude=[=](TdbEpoch t) {
        return CameraAttitude(offset(t.secondsPastJ2000-scene.epoch)*scene.nominalRotation);
    };
    synth::SyntheticImageOptions options;
    options.size={64,64}; options.subpixels=subpixels; options.exposureSeconds=.06; options.exposureSamples=5;
    options.fullPhaseElectrons=scene.budget; options.psfSigmaPixels=n.psfX; options.psfSigmaLinePixels=n.psfY;
    options.usePckEllipsoid=perturbed; options.albedoContrast=n.albedo; options.lommelSeeligerMix=n.reflectance;
    if(scene.controlledSun) options.sunDirectionCamera=offset(0)*(*scene.controlledSun);
    synth::OpticalDetectorNoise noise;
    noise.electronsPerDn=n.gain; noise.readNoiseElectrons=n.read; noise.biasDn=n.bias;
    noise.backgroundElectrons=50; noise.saturationDn=65535; noise.seed=scene.detectorSeed;
    Observation observed;
    observed.frame=synth::SyntheticImageGenerator(camera,attitude,options).createImage(
        "ENCELADUS",TdbEpoch(scene.epoch),noise,truth,"SIMULATED_CASSINI");
    observed.image.body="ENCELADUS"; observed.image.time=TdbEpoch(scene.epoch); observed.image.exposureSeconds=.06;
    observed.image.attitude=CameraAttitude(scene.nominalRotation);
    observed.image.dn=observed.frame.dn; observed.image.mask=observed.frame.mask;
    const CameraModel nominal(scene.intrinsics);
    observed.nominalPixel=nominal.project(scene.direction.vectorKm,observed.image.attitude);
    const fd::dynamics::SpiceStateProvider target("ENCELADUS");
    const auto pixel=[&](double t) {
        const auto d=core::computeApparentDirection(TdbEpoch(t),truth.stateAt(TdbEpoch(t)),target);
        return nominal.project(d.vectorKm,observed.image.attitude);
    };
    const auto first=pixel(scene.epoch-.03),last=pixel(scene.epoch+.03);
    observed.image.smearPixels={last.sample-first.sample,last.line-first.line};
    observed.calibration.intrinsics=scene.intrinsics; observed.calibration.saturationDn=65535;
    if(perturbed) {
        observed.image.attitudeCovariance=Eigen::Matrix3d::Identity()*std::pow(.25/scale,2);
        observed.calibration.principalPointCovariance=Eigen::Matrix2d::Identity()*.01;
        observed.calibration.relativeFocalSigma=.0005;
    }
    return observed;
}

BodyFitResult fit(const Scene& scene, const Observation& observed, bool ellipsoid, int order=24,
                  bool noiseless=false, int exposureSamples=1,
                  std::optional<BodyParameters> guess=std::nullopt) {
    BodyModel body;
    if(ellipsoid) body=makeSpiceBodyModel("ENCELADUS",scene.direction.lightTime.emissionEpoch,
        observed.image.attitude,scene.intrinsics);
    BodyFitOptions options; options.exposureSamples=exposureSamples;
    OpticalImage image=observed.image;
    if(noiseless) {
        const Nuisance n=ellipsoid ? scene.nuisance : Nuisance{};
        image.dn=(observed.frame.expectedSourceElectrons.array()+50)/n.gain+n.bias;
    }
    return BodyEstimator(BodyRenderer(body,order),observed.calibration,options).fit(image,guess);
}

void writeFit(std::ofstream& file, const Scene& s, const Observation& observed,
              const BodyFitResult& result, const std::string& model) {
    const bool converged=result.status==image::FitStatus::Converged && result.measurement.has_value();
    const Eigen::Matrix2d effective=converged ? result.measurement->covariance : Eigen::Matrix2d::Constant(missingValue);
    const auto center=result.parameters.center, actual=observed.frame.truth.geometricCenter;
    file<<s.experiment<<','<<model<<','<<s.id<<','<<s.draw<<','<<s.track<<','<<s.day<<','<<s.epoch<<','
        <<s.phaseDegrees<<','<<s.radiusPixels<<','<<s.budget<<','<<observed.frame.truth.phaseAngleRadians*180/std::numbers::pi<<','
        <<observed.frame.truth.radiusKm/observed.frame.truth.rangeKm*s.intrinsics.focalLength/.012<<','
        <<observed.frame.expectedSourceElectrons.sum()<<','<<int(result.status)<<','<<converged<<','<<result.passesFitChecks()<<','
        <<(converged ? center.sample-actual.sample : missingValue)<<','<<(converged ? center.line-actual.line : missingValue)<<','
        <<(converged ? center.sample-observed.nominalPixel.sample : missingValue)<<','<<(converged ? center.line-observed.nominalPixel.line : missingValue)<<','
        <<result.reducedChiSquared<<','<<result.formalCenterCovariance(0,0)<<','<<result.formalCenterCovariance(0,1)<<','
        <<result.formalCenterCovariance(1,1)<<','<<effective(0,0)<<','<<effective(0,1)<<','<<effective(1,1)<<','
        <<observed.frame.mask.cast<int>().sum()<<','<<observed.image.smearPixels.norm()<<'\n';
}

void writeRejection(std::ofstream& file,const Scene& scene,const Observation& observed,
                    const BodyFitResult& fit,const std::string& model) {
    BodyModel body;
    if(model=="ellipsoid")body=makeSpiceBodyModel("ENCELADUS",scene.direction.lightTime.emissionEpoch,
        observed.image.attitude,scene.intrinsics);
    const BodyRenderer renderer(body,8);
    RejectionContext context;context.calibration=&observed.calibration;context.renderer=&renderer;
    // No truth center, uncertainty, navigation innovation or truth-seeded alternatives.
    const auto decision=EstimatorRejection().evaluate(observed.image,fit,context);
    file<<scene.id<<','<<model<<','<<decision.accepted()<<',';
    for(std::size_t i=0;i<decision.reasons.size();++i)file<<(i?";":"")<<rejectionReasonName(decision.reasons[i]);
    for(const auto& value:{decision.peakSignalToNoise,decision.sourceSignalToNoise,decision.maskedFraction,
        decision.boundaryMarginPixels,decision.formalCenterSigmaPixels,decision.measurementCenterSigmaPixels})
        file<<','<<value.value_or(missingValue);
    file<<','<<int(decision.checks[static_cast<std::size_t>(RejectionCheck::Ambiguity)])<<','
        <<int(decision.checks[static_cast<std::size_t>(RejectionCheck::Navigation)])<<'\n';
}

int main(int argc,char** argv) {try {
    const fs::path root=DEEPNAV_SOURCE_DIR;
    fs::path directory=root/"results/cassini_synthetic_optical/validation";
    int draws=12,tracks=24,subpixels=8,selectedScene=-1;
    std::uint64_t seed=20261017;
    bool explicitOutput=false;
    for(int i=1;i<argc;++i) {
        const std::string key=argv[i];
        if(key=="--help") {std::cout<<"cassini_optical_validation_sweep [--output DIR] [--draws N] [--tracks N] [--subpixels N] [--scene N] [--seed N]\n";return 0;}
        if(++i==argc) throw std::invalid_argument("Missing option value.");
        if(key=="--output") {directory=fs::absolute(argv[i]);explicitOutput=true;}
        else if(key=="--draws") draws=std::stoi(argv[i]);
        else if(key=="--tracks") tracks=std::stoi(argv[i]);
        else if(key=="--subpixels") subpixels=std::stoi(argv[i]);
        else if(key=="--scene") selectedScene=std::stoi(argv[i]);
        else if(key=="--seed") seed=std::stoull(argv[i]);
        else throw std::invalid_argument("Unknown option.");
    }
    if(draws<1 || draws>1000 || tracks<1 || tracks>1000 || subpixels<4 || subpixels>16 || seed<1 || seed>1000000000)
        throw std::invalid_argument("Invalid experiment size or rendering resolution.");
    if(selectedScene>=27*draws+4*tracks || selectedScene< -1) throw std::invalid_argument("Replay scene outside experiment.");
    if(selectedScene>=0 && !explicitOutput) directory=root/("results/cassini_synthetic_optical/validation_replay_"+std::to_string(selectedScene));
    od::SpiceErrorModeGuard errors("RETURN","NONE");
    auto paths=fd::simulation::cassiniRadioKernels(root/"Kernels");
    for(const char* file:{"Cassini/fk/cas_v40.tf","Cassini/sclk/cas00172.tsc","Cassini/ck/04256_04292py_as_flown.bc",
        "Cassini/ck/04292_04320py_as_flown.bc","Cassini/ck/04281_04286ra.bc","Cassini/ck/04286_04291ra.bc",
        "Cassini/ck/04291_04296ra.bc","Cassini/ck/04296_04301ra.bc","Cassini/ck/04293_04298ra.bc"}) paths.push_back(root/"Kernels"/file);
    const fd::dynamics::SpiceKernelSet kernels(paths);
    const CassiniAttitude forceAttitude(paths);
    const Sat427Gravity gravity(root/"Kernels/Cassini/spk/200128RU_SCPSE_04251_04303.bsp.lbl");
    double recoilEpoch;str2et_c("2005-08-13T00:00:00",&recoilEpoch);
    double radii[3];SpiceInt count;bodvrd_c("SATURN","RADII",3,&count,radii);
    od::throwIfSpiceFailed("Sweep force setup");
    const fd::simulation::CassiniForceModel force(gravity,forceAttitude,fd::simulation::cassiniPerturbers(gravity),2,3220,12,1.42,
        {-2.725e-13,-2.725e-13,-5.232237415805135e-12},recoilEpoch,radii[0]);
    const fd::dynamics::SpiceStateProvider saturn("SATURN"),target("ENCELADUS");
    const fs::path truthPath=root/"results/cassini_two_way_wls/pass_schedule/truth/trajectory.csv";
    const fd::simulation::SavedSimulation truth(truthPath,force,saturn);
    bodvrd_c("ENCELADUS","RADII",3,&count,radii);od::throwIfSpiceFailed("Sweep target radii");
    const double meanRadius=(radii[0]+radii[1]+radii[2])/3;
    std::mt19937_64 random(seed);
    std::normal_distribution<double> normal;
    fs::create_directories(directory);
    auto results=output(directory/"observations.csv"),nuisances=output(directory/"nuisances.csv"),checks=output(directory/"resolution.csv");
    auto diagnostics=output(directory/"replay.csv");
    auto rejections=output(directory/"rejections.csv");
    rejections<<"scene,model,accepted,reasons,peak_snr,source_snr,masked_fraction,boundary_margin,formal_sigma,measurement_sigma,ambiguity_check,navigation_check\n";
    auto rejectionOptions=output(directory/"rejection_options.json");
    const RejectionOptions gates;
    rejectionOptions<<"{\"maximum_chi2\":"<<gates.maximumReducedChiSquared<<",\"minimum_peak_snr\":"<<gates.minimumPeakSignalToNoise
        <<",\"minimum_source_snr\":"<<gates.minimumSourceSignalToNoise<<",\"maximum_masked_fraction\":"<<gates.maximumMaskedFraction
        <<",\"minimum_valid_pixels\":"<<gates.minimumValidPixels<<",\"psf_support_sigmas\":"<<gates.psfSupportSigmas
        <<",\"minimum_boundary_margin\":"<<gates.minimumBoundaryMarginPixels<<",\"maximum_formal_sigma\":"<<gates.maximumFormalCenterSigmaPixels
        <<",\"maximum_measurement_sigma\":"<<gates.maximumMeasurementCenterSigmaPixels<<",\"ambiguity\":\"not supplied\",\"navigation\":\"not supplied\","
        <<"\"selection\":\"Image and nominal fit quantities only; no truth-based rejection or threshold tuning\"}\n";
    diagnostics<<"scene,model,initialization,status,accepted,extract_dx,extract_dy,radius_fit,phase_fit,chi2,common_window_chi2,truth_sample,truth_line\n";
    results<<"experiment,model,scene,draw,track,day,tdb,phase_requested,radius_requested,budget,phase_actual,radius_actual,source_electrons,status,converged,accepted,extract_dx,extract_dy,nav_dx,nav_dy,chi2,formal_ss,formal_sl,formal_ll,effective_ss,effective_sl,effective_ll,masked,nominal_smear_pixels\n";
    nuisances<<"scene,principal_x,principal_y,focal,k1,pointing_x,pointing_y,psf_x,psf_y,albedo,ls_mix,jitter,jitter_phase,gain,read,bias\n";
    checks<<"scene,model,phase,quad24_to_40_pixels,baseline_to_sub16_pixels,midpoint_to_three_sample_pixels,base_converged,quad_converged,fine_converged,integrated_converged\n";
    auto config=output(directory/"config.json");
    config<<"{\"draws_per_cell\":"<<draws<<",\"tracks\":"<<tracks<<",\"subpixels\":"<<subpixels
        <<",\"phases_degrees\":[30,75,120],\"radii_pixels\":[2,6,10],\"full_phase_electrons\":[20000,200000,2000000],"
        <<"\"selected_scene\":"<<selectedScene<<",\"grid_cells\":27,\"dates_utc\":[\"2004-10-10\",\"2004-10-11\",\"2004-10-12\",\"2004-10-13\"],"
        <<"\"seed\":"<<seed<<",\"fit_exposure_samples\":1,\"truth_exposure_samples\":5,\"exposure_seconds\":0.06,\"truth_csv\":"<<std::quoted(truthPath.string())<<",\"loaded_kernels\":[";
    for(std::size_t i=0;i<paths.size();++i) config<<(i?",":"")<<std::quoted(paths[i].string());
    config<<"],\"notes\":\"Synthetic controlled phase uses a Sun override; size uses focal-scale changes. All observer states use SavedSimulation. Sequences use physical SPK Sun and fixed focal scale. No EKF or uncertainty tuning.\"}\n";
    int sceneId=0,renderedScenes=0;
    const auto makeScene=[&](int day,double phase,double radius,double budget,int draw,int track,const std::string& experiment) {
        Scene s; s.id=sceneId++;s.day=day;s.draw=draw;s.track=track;s.experiment=experiment;
        s.detectorSeed=seed==20261017?800000+s.id:seed*1000003+s.id;
        const std::string utc="2004-10-"+std::to_string(10+day)+"T00:00:00";
        str2et_c(utc.c_str(),&s.epoch);od::throwIfSpiceFailed("Sweep epoch");
        if(s.epoch-.03<truth.startEpoch() || s.epoch+.03>truth.endEpoch()) throw std::runtime_error("Sweep outside truth coverage.");
        s.direction=core::computeApparentDirection(TdbEpoch(s.epoch),truth.stateAt(TdbEpoch(s.epoch)),target);
        const double scale=experiment=="grid" ? radius*s.direction.vectorKm.norm()/meanRadius : 2003.44/.012;
        s.intrinsics={scale*.012,{31.23,30.74},-Eigen::Matrix2d::Identity()/.012};
        const Eigen::Vector3d z=s.direction.unitDirection,x=z.unitOrthogonal(),y=z.cross(x);
        s.nominalRotation.row(0)=x.transpose();s.nominalRotation.row(1)=y.transpose();s.nominalRotation.row(2)=z.transpose();
        s.nominalRotation=Eigen::AngleAxisd(4/scale,Eigen::Vector3d::UnitX()).toRotationMatrix()*
            Eigen::AngleAxisd(-6/scale,Eigen::Vector3d::UnitY()).toRotationMatrix()*s.nominalRotation;
        s.phaseDegrees=phase;s.radiusPixels=radius;s.budget=budget;
        if(experiment=="grid") {
            const Eigen::Vector3d observer=-(s.nominalRotation*z).normalized();
            const Eigen::Vector3d a=(Eigen::Vector3d::UnitX()-observer.x()*observer).normalized(),b=observer.cross(a);
            const double azimuth=.4+day*.7,angle=phase*std::numbers::pi/180;
            s.controlledSun=std::cos(angle)*observer+std::sin(angle)*(std::cos(azimuth)*a+std::sin(azimuth)*b);
        }
        return s;
    };
    const auto run=[&](const Scene& s,bool includeControl) {
        if(selectedScene>=0 && selectedScene!=s.id) return;
        ++renderedScenes;
        const auto& n=s.nuisance;
        nuisances<<s.id<<','<<n.principal.x()<<','<<n.principal.y()<<','<<n.focal<<','<<n.k1<<','<<s.pointing.x()<<','<<s.pointing.y()<<','
            <<n.psfX<<','<<n.psfY<<','<<n.albedo<<','<<n.reflectance<<','<<n.jitter<<','<<n.jitterPhase<<','<<n.gain<<','<<n.read<<','<<n.bias<<'\n';
        std::optional<Observation> control;
        const auto record=[&](const Observation& observed,const BodyFitResult& fit,const std::string& model) {
            writeFit(results,s,observed,fit,model);writeRejection(rejections,s,observed,fit,model);
        };
        if(includeControl) {control=render(s,false,subpixels,truth);record(*control,fit(s,*control,false),"matched");}
        const auto realistic=render(s,true,subpixels,truth);
        record(realistic,fit(s,realistic,false),"sphere");
        const auto ellipsoid=fit(s,realistic,true);record(realistic,ellipsoid,"ellipsoid");
        // Replay only: truth-based initialization diagnoses difficult scenes;
        // it is never used in normal validation metrics or production fitting.
        if(selectedScene>=0) {
            for(const auto& model:{std::string("matched"),std::string("sphere"),std::string("ellipsoid")}) {
                if(model=="matched" && !control) continue;
                const auto& observed=model=="matched" ? *control : realistic;
                const double radius=observed.frame.truth.radiusKm/observed.frame.truth.rangeKm*s.intrinsics.focalLength/.012;
                const Eigen::Vector3d sun=observed.frame.truth.sunDirectionCamera;
                BodyParameters guess{observed.frame.truth.geometricCenter,radius,observed.frame.truth.phaseAngleRadians,
                    std::atan2(-sun.y(),-sun.x()),s.budget/(10*2*std::numbers::pi/3*radius*radius),105};
                for(bool seeded:{false,true}) {
                    const auto result=fit(s,observed,model=="ellipsoid",24,false,1,seeded?std::optional(guess):std::nullopt);
                    // Compare both solutions on exactly the same pixels, since
                    // an explicit initial center also changes the fitter's window.
                    const int x0=std::clamp(int(std::round(guess.center.sample))-18,0,27);
                    const int y0=std::clamp(int(std::round(guess.center.line))-18,0,27);
                    BodyModel comparisonBody;
                    if(model=="ellipsoid") comparisonBody=makeSpiceBodyModel("ENCELADUS",s.direction.lightTime.emissionEpoch,
                        observed.image.attitude,s.intrinsics);
                    const auto predicted=BodyRenderer(comparisonBody,24).render({37,37},{double(x0),double(y0)},result.parameters,.8).dn;
                    const image::Image pixels=observed.image.dn.block(y0,x0,37,37);
                    const image::Image variance=(pixels.array()-100).max(0)/10+std::pow(3./10,2)+1./12;
                    const double commonChi=((pixels-predicted).array().square()/variance.array()).sum()/(37*37-7);
                    diagnostics<<s.id<<','<<model<<','<<(seeded?"truth_seed_diagnostic":"automatic")<<','<<int(result.status)<<','<<result.passesFitChecks()<<','
                        <<result.parameters.center.sample-observed.frame.truth.geometricCenter.sample<<','
                        <<result.parameters.center.line-observed.frame.truth.geometricCenter.line<<','<<result.parameters.radiusPixels<<','
                        <<result.parameters.phaseRadians*180/std::numbers::pi<<','<<result.reducedChiSquared<<','<<commonChi<<','
                        <<observed.frame.truth.geometricCenter.sample<<','<<observed.frame.truth.geometricCenter.line<<'\n';
                }
                auto pixels=output(directory/(model+"_replay.csv"));
                pixels<<"sample,line,dn,source_electrons\n";
                for(int y=0;y<64;++y)for(int x=0;x<64;++x) pixels<<x<<','<<y<<','<<observed.image.dn(y,x)<<','<<observed.frame.expectedSourceElectrons(y,x)<<'\n';
            }
        }
        // Export nine representative images using identical scales in the Python montage.
        if(s.experiment=="grid" && s.draw==0 && s.radiusPixels==6) {
            auto pixels=output(directory/("scene_"+std::to_string(s.id)+".csv"));
            pixels<<"sample,line,dn,source_electrons,model_dn,residual_sigma\n";
            for(int y=0;y<64;++y)for(int x=0;x<64;++x) {
                const int xx=x-int(ellipsoid.origin.sample),yy=y-int(ellipsoid.origin.line);
                const bool inside=xx>=0 && yy>=0 && xx<ellipsoid.modelDn.cols() && yy<ellipsoid.modelDn.rows();
                pixels<<x<<','<<y<<','<<realistic.image.dn(y,x)<<','<<realistic.frame.expectedSourceElectrons(y,x)<<','
                    <<(inside?ellipsoid.modelDn(yy,xx):missingValue)<<','<<(inside?ellipsoid.residualDn(yy,xx)/std::sqrt(ellipsoid.varianceDn2(yy,xx)):missingValue)<<'\n';
            }
        }
        // Independent noiseless checks distinguish image/fit numerical error from detector noise.
        if(s.experiment=="grid" && s.draw==0 && s.radiusPixels==6 && s.budget==200000) {
            for(bool perturbed:{false,true}) {
                const auto& coarse=perturbed?realistic:*control;
                const auto fine=render(s,perturbed,16,truth);
                const auto base=fit(s,coarse,perturbed,24,true),quad=fit(s,coarse,perturbed,40,true),refined=fit(s,fine,perturbed,40,true);
                const auto integrated=fit(s,coarse,perturbed,24,true,3);
                const auto distance=[](const BodyFitResult& a,const BodyFitResult& b) {
                    if(!a.measurement || !b.measurement)return missingValue;
                    return std::hypot(a.parameters.center.sample-b.parameters.center.sample,a.parameters.center.line-b.parameters.center.line);
                };
                checks<<s.id<<','<<(perturbed?"ellipsoid":"matched")<<','<<s.phaseDegrees<<','<<distance(base,quad)<<','<<distance(quad,refined)<<','<<distance(base,integrated)<<','
                    <<bool(base.measurement)<<','<<bool(quad.measurement)<<','<<bool(refined.measurement)<<','<<bool(integrated.measurement)<<'\n';
            }
        }
        results.flush();rejections.flush();nuisances.flush();checks.flush();
    };
    for(double phase:{30.,75.,120.})for(double radius:{2.,6.,10.})for(double budget:{20000.,200000.,2000000.}) {
        for(int draw=0;draw<draws;++draw) {
            auto s=makeScene(draw%4,phase,radius,budget,draw,-1,"grid");
            s.nuisance=drawNuisance(random,.25);s.pointing=s.nuisance.alignment;
            run(s,true);
        }
        if(selectedScene<0) std::cout<<"Grid phase="<<phase<<" radius="<<radius<<" budget="<<budget<<": "<<draws<<" independent draws complete\n"<<std::flush;
    }
    for(int track=0;track<tracks;++track) {
        const auto shared=drawNuisance(random,.20);
        for(int day=0;day<4;++day) {
            auto s=makeScene(day,-1,-1,200000,day,track,"sequence");s.nuisance=shared;
            s.nuisance.jitter=std::abs(.08*normal(random));s.nuisance.jitterPhase=double(day);
            s.pointing=shared.alignment+Eigen::Vector2d(.15*normal(random),.15*normal(random));
            run(s,false);
        }
        if(selectedScene<0) std::cout<<"Shared-bias sequence "<<track+1<<'/'<<tracks<<" complete\n"<<std::flush;
    }
    std::cout<<"Saved "<<renderedScenes<<" scenes to "<<directory<<'\n';
} catch(const std::exception& error) {std::cerr<<"ERROR: "<<error.what()<<'\n';return 1;} }
