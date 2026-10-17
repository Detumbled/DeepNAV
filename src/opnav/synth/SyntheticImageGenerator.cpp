#include "opnav/synth/SyntheticImageGenerator.hpp"
#include "dynamics/SpiceStateProvider.hpp"
#include "opnav/core/ApparentDirection.hpp"
#include "utils/CSPICE/SpiceError.hpp"
#include <Eigen/LU>
#include <algorithm>
#include <cmath>
#include <numbers>
#include <random>
#include <stdexcept>
#include <utility>

namespace fd::opnav::synth {
namespace {
void require(bool condition, const char* message) {
    if (!condition) throw std::invalid_argument(message);
}
bool nonnegative(double x) { return std::isfinite(x) && x >= 0; }
std::optional<PixelCoordinates> centerOfLight(const image::Image& source) {
    const double sum = source.sum();
    if (!(sum > 0)) return std::nullopt;
    double x = 0, y = 0;
    for (int l = 0; l < source.rows(); ++l) for (int s = 0; s < source.cols(); ++s) {
        x += s * source(l,s); y += l * source(l,s);
    }
    return PixelCoordinates{x/sum, y/sum};
}
// Unit-sum discrete Gaussian on the oversampled grid. Zero padding preserves
// loss outside the image; it does not renormalize a clipped target's flux.
image::Image blur(const image::Image& source,double sigmaX,double sigmaY) {
    const auto kernel=[](double sigma) -> Eigen::VectorXd {
        if(sigma==0)return Eigen::VectorXd::Ones(1);
        const int radius=int(std::ceil(5*sigma));Eigen::VectorXd values(2*radius+1);
        for(int k=-radius;k<=radius;++k)values[k+radius]=std::exp(-.5*k*k/(sigma*sigma));
        return values/values.sum();
    };
    const Eigen::VectorXd kx=kernel(sigmaX),ky=kernel(sigmaY);
    const int rx=kx.size()/2,ry=ky.size()/2;
    image::Image temporary=image::Image::Zero(source.rows(),source.cols()),result=temporary;
    for(int y=0;y<source.rows();++y)for(int k=-rx;k<=rx;++k) {
        const int begin=std::max(0,-k),end=std::min(int(source.cols()),int(source.cols())-k);
        if(end>begin)temporary.row(y).segment(begin,end-begin)+=kx[k+rx]*source.row(y).segment(begin+k,end-begin);
    }
    for(int k=-ry;k<=ry;++k) {
        const int begin=std::max(0,-k),end=std::min(int(source.rows()),int(source.rows())-k);
        if(end>begin)result.middleRows(begin,end-begin)+=ky[k+ry]*temporary.middleRows(begin+k,end-begin);
    }
    return result;
}
}
struct SyntheticImageGenerator::Geometry {
    OpticalFrameTruth truth;
    Eigen::Vector3d centerCameraKm, sunCamera;
};

SyntheticImageGenerator::SyntheticImageGenerator(CameraModel camera, AttitudeProvider attitude,
                                                 SyntheticImageOptions options)
    : camera_(std::move(camera)), attitude_(std::move(attitude)), options_(std::move(options)) {
    require(bool(attitude_), "An independent truth attitude provider is required.");

    const auto& k = camera_.getIntrinsics();
    require(std::isfinite(k.focalLength) && k.focalLength > 0 &&
            std::isfinite(k.principalPoint.sample) && std::isfinite(k.principalPoint.line) &&
            k.focalPlaneToPixelMatrix_PixelPermm.allFinite() &&
            std::abs(k.focalPlaneToPixelMatrix_PixelPermm.determinant()) > 1e-15, "Invalid pinhole calibration.");
    require(options_.size.width > 0 && options_.size.height > 0 && options_.size.width <= 2048 && options_.size.height <= 2048 &&
            options_.subpixels >= 1 && options_.subpixels <= 16 &&
            std::int64_t(options_.size.width)*options_.size.height*options_.subpixels*options_.subpixels <= 16777216,
            "Invalid or excessive synthetic image grid.");
    require(nonnegative(options_.radiusKm) && nonnegative(options_.fullPhaseElectrons) &&
            options_.fullPhaseElectrons <= 1e12 &&
            nonnegative(options_.psfSigmaPixels) && options_.psfSigmaPixels <= 10 &&
            nonnegative(options_.exposureSeconds) && options_.exposureSamples >= 1 && options_.exposureSamples <= 101,
            "Invalid sphere, flux, PSF, or exposure options.");
    require(nonnegative(options_.psfSigmaLinePixels.value_or(options_.psfSigmaPixels)) && options_.psfSigmaLinePixels.value_or(0)<=10 &&
            nonnegative(options_.albedoContrast) && options_.albedoContrast<1 && nonnegative(options_.lommelSeeligerMix) && options_.lommelSeeligerMix<=1 &&
            nonnegative(options_.referenceRangeKm) && options_.referenceExposureSeconds>0 && std::isfinite(options_.referenceExposureSeconds), "Invalid realism options.");
    if(options_.bodyAxesKm)require(options_.bodyAxesKm->allFinite() && (options_.bodyAxesKm->array()>0).all(),"Invalid body axes.");
    const int padding = int(std::ceil(5*std::max(options_.psfSigmaPixels,options_.psfSigmaLinePixels.value_or(0))))+1;
    require(std::int64_t(options_.size.width+2*padding)*(options_.size.height+2*padding)*
            options_.subpixels*options_.subpixels <= 16777216, "Padded image grid exceeds memory budget.");
    if (options_.sunDirectionCamera) require(options_.sunDirectionCamera->allFinite() &&
        options_.sunDirectionCamera->norm() > 0, "Invalid synthetic illumination direction.");
}

SyntheticImageGenerator::Geometry SyntheticImageGenerator::geometry(const std::string& body,
    TdbEpoch epoch, const dynamics::StateProvider& observerProvider, const std::string& observer) const {
    require(!body.empty() && !observer.empty() && body != observer && std::isfinite(epoch.secondsPastJ2000),
            "Invalid target, observer, or TDB epoch.");
    require(observerProvider.referenceFrame()=="J2000" && observerProvider.referenceOrigin()=="SSB",
            "Observer truth must provide geometric J2000/SSB position and velocity.");
    const dynamics::SpiceStateProvider targetProvider(body), sunProvider("SUN");
    const auto observerState = observerProvider.stateAt(epoch);
    const auto direction = core::computeApparentDirection(epoch, observerState, targetProvider);
    const auto attitude = attitude_(epoch);
    const double range = direction.lightTime.geometricLineOfSightKm.norm();
    Eigen::Vector3d axes=Eigen::Vector3d::Zero();
    double radius = options_.radiusKm;
    if (radius == 0 || options_.usePckEllipsoid) {
        od::SpiceErrorModeGuard guard;
        SpiceDouble radii[3]; SpiceInt count;
        bodvrd_c(body.c_str(), "RADII", 3, &count, radii);
        od::throwIfSpiceFailed("Synthetic target PCK radii");
        require(count == 3, "Target needs three PCK radii or an explicit sphere radius.");
        axes=Eigen::Map<Eigen::Vector3d>(radii);
        if(radius==0)radius=axes.mean();
    }
    require(std::isfinite(radius) && radius > 0 && radius/range < 0.01,
            "Renderer requires a distant sphere (radius/range < 0.01).");
    if(!options_.usePckEllipsoid)axes.setConstant(radius);
    if(options_.bodyAxesKm)axes=*options_.bodyAxesKm;
    require(axes.maxCoeff()/range<0.01,"Renderer requires a distant ellipsoid.");
    Eigen::Matrix3d bodyToCamera=Eigen::Matrix3d::Identity();
    if(options_.usePckEllipsoid || options_.bodyAxesKm || options_.albedoContrast>0) {
        od::SpiceErrorModeGuard guard;SpiceInt frameId;SpiceChar name[128];SpiceBoolean found;
        cnmfrm_c(body.c_str(),sizeof(name),&frameId,name,&found);od::throwIfSpiceFailed("Target body-fixed frame");
        require(found,"No body-fixed frame for target.");SpiceDouble rotation[3][3];
        pxform_c(name,"J2000",direction.lightTime.emissionEpoch.secondsPastJ2000,rotation);od::throwIfSpiceFailed("Body orientation at emission");
        Eigen::Matrix3d inertial;for(int i=0;i<3;++i)for(int j=0;j<3;++j)inertial(i,j)=rotation[i][j];
        bodyToCamera=attitude.matrix()*inertial;
    }
    const Eigen::Vector3d center = attitude.toCamera(direction.unitDirection)*range;
    const auto pixel = camera_.project(center);
    Eigen::Vector3d sun;
    if (options_.sunDirectionCamera) sun = options_.sunDirectionCamera->normalized();
    else {
        // Geometric illumination at target emission epoch; solar light time and
        // relativistic surface effects are outside this simple brightness model.
        const auto solar = sunProvider.stateAt(direction.lightTime.emissionEpoch);
        sun = attitude.toCamera(solar.positionKm-direction.lightTime.targetStateAtEmission.positionKm).normalized();
    }
    require(sun.allFinite() && sun.norm() > 0, "Cannot determine a finite illumination direction.");
    OpticalFrameTruth truth;
    truth.target = body; truth.observer = observer; truth.epoch = epoch;
    truth.geometricCenter = pixel; truth.rangeKm = range; truth.radiusKm = radius;
    truth.lightTimeSeconds = direction.lightTime.lightTimeSeconds;
    truth.phaseAngleRadians = std::acos(std::clamp(sun.dot(-center.normalized()),-1.0,1.0));
    truth.observerPositionKm = observerState.positionKm; truth.observerVelocityKmPerSecond = observerState.velocityKmPerSec;
    truth.inertialToCamera = attitude.matrix();truth.bodyAxesKm=axes;truth.bodyToCamera=bodyToCamera;truth.sunDirectionCamera=sun;
    return {std::move(truth), center, sun};
}

SyntheticOpticalMeasurement SyntheticImageGenerator::createMeasurement(const std::string& body,
    TdbEpoch tdb, OpticalPixelNoise noise, const std::string& observer) const {
    return createMeasurement(body,tdb,noise,dynamics::SpiceStateProvider(observer),observer);
}

SyntheticOpticalMeasurement SyntheticImageGenerator::createMeasurement(const std::string& body,
    TdbEpoch tdb, OpticalPixelNoise noise, const dynamics::StateProvider& observer,
    const std::string& observerLabel) const {
    require(nonnegative(noise.sigmaPixels) && noise.sigmaPixels <= 1e6, "Invalid pixel noise sigma.");
    auto g = geometry(body,tdb,observer,observerLabel);
    auto measured = g.truth.geometricCenter;
    require(measured.sample >= -0.5 && measured.sample < options_.size.width-0.5 &&
            measured.line >= -0.5 && measured.line < options_.size.height-0.5,
            "Geometric center is outside the synthetic detector.");
    if (noise.sigmaPixels > 0) {
        std::mt19937_64 engine(noise.seed);
        std::normal_distribution<double> distribution(0,noise.sigmaPixels);
        measured.sample += distribution(engine); measured.line += distribution(engine);
    }
    g.truth.seed = noise.seed;
    return {{measured, Eigen::Matrix2d::Identity()*noise.sigmaPixels*noise.sigmaPixels},std::move(g.truth)};
}

SyntheticOpticalFrame SyntheticImageGenerator::createImage(const std::string& body, TdbEpoch tdb,
    OpticalDetectorNoise noise, const std::string& observer) const {
    return createImage(body,tdb,noise,dynamics::SpiceStateProvider(observer),observer);
}

SyntheticOpticalFrame SyntheticImageGenerator::createImage(const std::string& body, TdbEpoch tdb,
    OpticalDetectorNoise noise, const dynamics::StateProvider& observer, const std::string& observerLabel) const {
    require(std::isfinite(noise.electronsPerDn) && noise.electronsPerDn >= 1e-12 &&
            nonnegative(noise.backgroundElectrons) && noise.backgroundElectrons <= 1e12 &&
            nonnegative(noise.readNoiseElectrons) && noise.readNoiseElectrons <= 1e6 &&
            nonnegative(noise.biasDn) && std::isfinite(noise.saturationDn) && noise.saturationDn > noise.biasDn,
            "Invalid detector configuration.");
    auto mid = geometry(body,tdb,observer,observerLabel);
    const auto& intrinsics = camera_.getIntrinsics();
    const Eigen::Matrix2d focal = intrinsics.focalLength*intrinsics.focalPlaneToPixelMatrix_PixelPermm;

    const Eigen::Vector3d c = mid.centerCameraKm;
    // Local projected-area Jacobian. Full-phase Lambert disk integral is 2*pi*R^2/3.
    Eigen::Matrix<double,2,3> perspective;
    perspective << 1/c.z(),0,-c.x()/(c.z()*c.z()), 0,1/c.z(),-c.y()/(c.z()*c.z());
    const Eigen::Vector3d z = c.normalized(), x = z.unitOrthogonal(), y = z.cross(x);
    Eigen::Matrix<double,3,2> tangent; tangent.col(0)=x; tangent.col(1)=y;
    const double area = std::abs((focal*perspective*tangent).determinant());
    const double budget=options_.fullPhaseElectrons*(options_.referenceRangeKm>0 ?
        std::pow(options_.referenceRangeKm/mid.truth.rangeKm,2)*options_.exposureSeconds/options_.referenceExposureSeconds : 1);
    const double amplitude = budget/(2*std::numbers::pi/3*mid.truth.radiusKm*mid.truth.radiusKm*area);
    require(std::isfinite(amplitude), "Invalid source radiance normalization.");
    const int sub = options_.subpixels, width = options_.size.width, height = options_.size.height;
    const int padding = int(std::ceil(5*std::max(options_.psfSigmaPixels,options_.psfSigmaLinePixels.value_or(0))))+1;
    // Render beyond the detector so off-frame light can still blur into it.
    image::Image high = image::Image::Zero((height+2*padding)*sub,(width+2*padding)*sub);
    Eigen::Matrix<double,Eigen::Dynamic,3,Eigen::RowMajor> rays(high.size(),3);
    for(int l=0;l<high.rows();++l)for(int s=0;s<high.cols();++s)
        rays.row(l*high.cols()+s)=camera_.unproject({(s+.5)/sub-padding-.5,(l+.5)/sub-padding-.5}).transpose();
    const int samples = options_.exposureSeconds == 0 ? 1 : options_.exposureSamples;
    for (int sample = 0; sample < samples; ++sample) {
        const double offset = options_.exposureSeconds*((sample+0.5)/samples-0.5);
        const auto g = geometry(body,TdbEpoch(tdb.secondsPastJ2000+offset),observer,observerLabel);
        const double r2 = g.truth.radiusKm*g.truth.radiusKm;
        const bool ellipsoid=options_.usePckEllipsoid || options_.bodyAxesKm.has_value();
        const Eigen::Matrix3d q=g.truth.bodyToCamera*g.truth.bodyAxesKm.cwiseInverse().cwiseAbs2().asDiagonal()*g.truth.bodyToCamera.transpose();
        for (int l = 0; l < high.rows(); ++l) for (int s = 0; s < high.cols(); ++s) {
            const Eigen::Vector3d ray=rays.row(l*high.cols()+s).transpose();
            double distance;
            if(ellipsoid) {
                const double a=ray.dot(q*ray),closest=ray.dot(q*g.centerCameraKm)/a;
                const Eigen::Vector3d delta=closest*ray-g.centerCameraKm;
                const double discriminant=1-delta.dot(q*delta);
                if(closest<=0 || discriminant<=0)continue;
                distance=closest-std::sqrt(discriminant/a);
            } else {
                const double along=ray.dot(g.centerCameraKm),discriminant=r2-g.centerCameraKm.cross(ray).squaredNorm();
                if(along<=0 || discriminant<=0)continue;
                distance=along-std::sqrt(discriminant);
            }
            const Eigen::Vector3d surface=distance*ray-g.centerCameraKm;
            const Eigen::Vector3d normal=ellipsoid ? (q*surface).normalized() : surface/g.truth.radiusKm;
            const double incidence=std::max(0.,normal.dot(g.sunCamera)),emission=std::max(0.,-normal.dot(ray));
            double brightness=(1-options_.lommelSeeligerMix)*incidence+options_.lommelSeeligerMix*incidence/std::max(1e-12,incidence+emission);
            if(options_.albedoContrast>0) {
                const Eigen::Vector3d body=(g.truth.bodyToCamera.transpose()*surface).cwiseQuotient(g.truth.bodyAxesKm).normalized();
                brightness*=1+options_.albedoContrast*std::cos(3*std::atan2(body.y(),body.x())+.7)*std::sqrt(std::max(0.,1-body.z()*body.z()));
            }
            high(l,s)+=amplitude*brightness/(sub*sub*samples);
        }
    }
    const auto convolved = blur(high,options_.psfSigmaPixels*sub,options_.psfSigmaLinePixels.value_or(options_.psfSigmaPixels)*sub);
    image::Image source(height,width);
    for (int l = 0; l < height; ++l) for (int s = 0; s < width; ++s)
        source(l,s) = convolved.block((l+padding)*sub,(s+padding)*sub,sub,sub).sum();
    mid.truth.brightnessCenter = centerOfLight(source); mid.truth.seed = noise.seed;
    SyntheticOpticalFrame result;
    result.expectedSourceElectrons = source; result.truth = std::move(mid.truth);
    auto detector=sampleDetector(source,noise);
    result.dn=std::move(detector.dn);result.varianceDn2=std::move(detector.varianceDn2);result.mask=std::move(detector.mask);
    return result;
}
} // namespace fd::opnav::synth
