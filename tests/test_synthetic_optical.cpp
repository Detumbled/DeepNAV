#include "opnav/synth/SyntheticImageGenerator.hpp"
#include "dynamics/SpiceKernelSet.hpp"
#include "dynamics/SpiceStateProvider.hpp"
#include "opnav/core/ApparentDirection.hpp"
#include "utils/CSPICE/SpiceError.hpp"
#include <cmath>
#include <filesystem>
#include <iostream>
#include <limits>
#include <numbers>
#include <stdexcept>
using namespace fd::opnav;
namespace {
void require(bool condition,const char* message) { if(!condition) throw std::runtime_error(message); }
template<class F> void rejects(F f) {
    bool rejected=false; try { f(); } catch(const std::exception&) { rejected=true; }
    require(rejected,"Invalid request was accepted.");
}
class OffsetObserver final : public fd::dynamics::StateProvider {
    const fd::dynamics::StateProvider& base_;
public:
    Eigen::Vector3d offset;
    bool invalidOrigin{false};
    OffsetObserver(const fd::dynamics::StateProvider& base, Eigen::Vector3d value):base_(base),offset(value){}
    fd::dynamics::CartesianState stateAt(TdbEpoch t) const override {
        auto state=base_.stateAt(t);state.positionKm+=offset;return state;
    }
    std::string_view referenceFrame() const noexcept override {return "J2000";}
    std::string_view referenceOrigin() const noexcept override {return invalidOrigin ? "SATURN" : "SSB";}
};
}
int main() {
    try {
        const std::filesystem::path root=DEEPNAV_SOURCE_DIR;
        const auto meta=root/"benchmarks/cassini_synthetic_optical/kernels.tm";
        if(!std::filesystem::exists(root/"Kernels/Cassini/spk/200128RU_SCPSE_04251_04303.bsp")) {
            std::cout << "Local Cassini kernels unavailable.\n"; return 77;
        }
        std::filesystem::current_path(root);
        od::SpiceErrorModeGuard guard("RETURN","NONE");
        const fd::dynamics::SpiceKernelSet kernels({meta});
        double et; str2et_c("2004-10-10T18:13:06.046",&et); od::throwIfSpiceFailed("Test epoch");
        const TdbEpoch epoch(et);
        const fd::dynamics::SpiceStateProvider observer("CASSINI"),target("ENCELADUS");
        const auto dir=core::computeApparentDirection(epoch,observer.stateAt(epoch),target);
        const Eigen::Vector3d z=dir.unitDirection,x=z.unitOrthogonal(),y=z.cross(x);
        Eigen::Matrix3d rotation; rotation.row(0)=x.transpose(); rotation.row(1)=y.transpose(); rotation.row(2)=z.transpose();
        const CameraAttitude attitude(rotation);
        const auto provider=[attitude](TdbEpoch){return attitude;};
        const CameraModel camera({2003.44,{31.5,31.5},-Eigen::Matrix2d::Identity()/0.012});
        synth::SyntheticImageOptions options;
        options.size={64,64}; options.subpixels=8; options.exposureSeconds=0; options.psfSigmaPixels=0;
        options.sunDirectionCamera=Eigen::Vector3d(0,0,-1);
        synth::OpticalDetectorNoise noiseless;
        noiseless.shotNoise=false; noiseless.quantize=false; noiseless.readNoiseElectrons=0;
        noiseless.backgroundElectrons=0; noiseless.saturationDn=1e12;
        synth::SyntheticImageGenerator full(camera,provider,options);
        const auto frame=full.create("ENCELADUS",epoch,noiseless,"CASSINI");
        const auto exact=full.createMeasurement("ENCELADUS",epoch,{0,1},"CASSINI");
        require(std::hypot(exact.measurement.center.sample-31.5,exact.measurement.center.line-31.5)<1e-8 &&
                exact.measurement.covariance.norm()==0,"Perfect direct observation is not exact.");
        OffsetObserver simulated(observer,100*x);
        const auto injected=full.createMeasurement("ENCELADUS",epoch,{0,1},simulated);
        const auto injectedFrame=full.create("ENCELADUS",epoch,noiseless,simulated);
        require((injectedFrame.truth.observerPositionKm-simulated.stateAt(epoch).positionKm).norm()<1e-10 &&
                (injectedFrame.truth.observerVelocityKmPerSecond-simulated.stateAt(epoch).velocityKmPerSec).norm()<1e-12,
                "Injected simulated observer state was replaced by SPK.");
        require(std::hypot(injected.measurement.center.sample-exact.measurement.center.sample,
                           injected.measurement.center.line-exact.measurement.center.line)>2,
                "Injected observer displacement did not change the optical observation.");
        simulated.invalidOrigin=true;
        rejects([&]{ (void)full.createMeasurement("ENCELADUS",epoch,{0,1},simulated); });
        require(std::abs(frame.expectedSourceElectrons.sum()/options.fullPhaseElectrons-1)<0.003,"Full-phase flux normalization failed.");
        require(frame.truth.brightnessCenter.has_value(),"Missing noiseless brightness center.");
        require(std::hypot(frame.truth.brightnessCenter->sample-31.5,frame.truth.brightnessCenter->line-31.5)<1e-6,"Symmetric source centroid is biased.");
        options.psfSigmaPixels=1.2;
        const auto blurred=synth::SyntheticImageGenerator(camera,provider,options).create("ENCELADUS",epoch,noiseless,"CASSINI");
        require(std::abs(blurred.expectedSourceElectrons.sum()/frame.expectedSourceElectrons.sum()-1)<1e-10,"PSF does not conserve enclosed flux.");
        // Exposure integration must preserve source motion instead of re-aiming each sample.
        options.exposureSeconds=30; options.exposureSamples=31; options.subpixels=4;
        const auto smeared=synth::SyntheticImageGenerator(camera,provider,options).create("ENCELADUS",epoch,noiseless,"CASSINI");
        auto lineVariance=[](const synth::SyntheticOpticalFrame& f) {
            double variance=0;
            for(int l=0;l<f.expectedSourceElectrons.rows();++l)
                variance+=std::pow(l-f.truth.brightnessCenter->line,2)*f.expectedSourceElectrons.row(l).sum();
            return variance/f.expectedSourceElectrons.sum();
        };
        require(lineVariance(smeared)>lineVariance(blurred)+0.4,"Exposure motion was not integrated.");
        options.exposureSeconds=0; options.subpixels=8;
        options.psfSigmaPixels=0; options.sunDirectionCamera=Eigen::Vector3d(1,0,0);
        const auto half=synth::SyntheticImageGenerator(camera,provider,options).create("ENCELADUS",epoch,noiseless,"CASSINI");
        require(std::abs(half.expectedSourceElectrons.sum()/frame.expectedSourceElectrons.sum()-1/std::numbers::pi)<0.004,"Lambert quarter-phase flux disagrees with analytic phase law.");
        require(half.truth.brightnessCenter->sample<31 && std::abs(half.truth.brightnessCenter->line-31.5)<0.01,"Illumination photocenter direction is incorrect.");
        options.sunDirectionCamera=Eigen::Vector3d(0,0,-1);
        const CameraModel edgeCamera({2003.44,{0.5,31.5},-Eigen::Matrix2d::Identity()/0.012});
        const auto clipped=synth::SyntheticImageGenerator(edgeCamera,provider,options).create("ENCELADUS",epoch,noiseless,"CASSINI");
        require(clipped.expectedSourceElectrons.sum()<0.8*frame.expectedSourceElectrons.sum(),"Clipped target flux was incorrectly renormalized.");
        auto noisy=noiseless; noisy.shotNoise=true; noisy.readNoiseElectrons=3; noisy.backgroundElectrons=50; noisy.seed=123;
        const auto a=full.create("ENCELADUS",epoch,noisy,"CASSINI"),b=full.create("ENCELADUS",epoch,noisy,"CASSINI");
        require((a.dn-b.dn).norm()==0,"Seeded images are not reproducible.");
        noisy.saturationDn=30;
        const auto saturated=full.create("ENCELADUS",epoch,noisy,"CASSINI");
        require(saturated.mask.cast<int>().sum()>0 && saturated.dn.maxCoeff()<=30,"Saturation/mask failed.");
        // Empty-source frames give independent detector statistics, without model-fit assumptions.
        options.size={16,16}; options.subpixels=1; options.fullPhaseElectrons=0;
        const synth::SyntheticImageGenerator blank(camera,provider,options);
        noisy.backgroundElectrons=100; noisy.electronsPerDn=2; noisy.biasDn=100; noisy.readNoiseElectrons=5; noisy.saturationDn=4095;
        double sum=0,squares=0; int n=0;
        for(int i=0;i<100;++i) {
            noisy.seed=400+i;
            const auto f=blank.create("ENCELADUS",epoch,noisy,"CASSINI");
            require(std::abs(f.varianceDn2(0,0)-31.25)<1e-12,"Detector variance is wrong.");
            for(int l=0;l<16;++l)for(int s=0;s<16;++s) { const double e=f.dn(l,s)-150; sum+=e; squares+=e*e; ++n; }
        }
        require(std::abs(sum/n)<0.2 && std::abs(squares/n/31.25-1)<0.05,"Poisson/read noise statistics failed.");
        double ps=0,pl=0,p2=0;
        for(int i=0;i<2000;++i) {
            const auto measurement=full.createMeasurement("ENCELADUS",epoch,{0.2,std::uint64_t(10000+i)},"CASSINI");
            const double dx=measurement.measurement.center.sample-measurement.truth.geometricCenter.sample;
            const double dy=measurement.measurement.center.line-measurement.truth.geometricCenter.line;
            ps+=dx; pl+=dy; p2+=dx*dx+dy*dy;
            require(std::abs(measurement.measurement.covariance(0,0)-0.04)<1e-12,"Pixel covariance is wrong.");
        }
        require(std::abs(ps/2000)<0.02 && std::abs(pl/2000)<0.02 && std::abs(p2/2000/0.08-1)<0.1,"Direct pixel noise/covariance inconsistent.");
        rejects([&]{ (void)full.createMeasurement("ENCELADUS",epoch,{-1,1},"CASSINI"); });
        rejects([&]{ (void)full.createMeasurement("ENCELADUS",TdbEpoch(std::numeric_limits<double>::quiet_NaN()),{},"CASSINI"); });
        rejects([&]{ (void)full.createMeasurement("ENCELADUS",epoch,{},"ENCELADUS"); });
        rejects([&]{ (void)full.createMeasurement("NOT_A_SPICE_BODY",epoch,{},"CASSINI"); });
        rejects([&]{ (void)full.createMeasurement("ENCELADUS",TdbEpoch(-1e10),{},"CASSINI"); });
        options.subpixels=0;
        rejects([&]{ synth::SyntheticImageGenerator bad(camera,provider,options); });
        std::cout << "PASS: symmetry, Lambert phase law, flux/PSF, clipping, seed reproducibility, saturation, detector and pixel statistics, invalid inputs and kernel coverage.\n";
    } catch(const std::exception& e) {std::cerr << "FAIL: " << e.what() << '\n';return 1;}
}
