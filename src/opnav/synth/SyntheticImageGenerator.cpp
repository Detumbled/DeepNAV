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
image::Image blur(const image::Image& source, double sigma) {
    if (sigma == 0) return source;
    const int radius = int(std::ceil(5*sigma));
    Eigen::VectorXd kernel(2*radius+1);
    for (int i = -radius; i <= radius; ++i) kernel[i+radius] = std::exp(-0.5*i*i/(sigma*sigma));
    kernel /= kernel.sum();
    image::Image temporary = image::Image::Zero(source.rows(),source.cols());
    image::Image result = temporary;
    for (int y = 0; y < source.rows(); ++y) for (int k = -radius; k <= radius; ++k) {
        const int begin = std::max(0,-k), end = std::min(int(source.cols()),int(source.cols())-k);
        if (end > begin) temporary.row(y).segment(begin,end-begin) += kernel[k+radius]*source.row(y).segment(begin+k,end-begin);
    }
    for (int k = -radius; k <= radius; ++k) {
        const int begin = std::max(0,-k), end = std::min(int(source.rows()),int(source.rows())-k);
        if (end > begin) result.middleRows(begin,end-begin) += kernel[k+radius]*temporary.middleRows(begin+k,end-begin);
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
    require(dynamic_cast<const NoDistortion*>(&camera_.getDistortionModel()), "Synthetic rendering currently requires NoDistortion.");
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
    const int padding = int(std::ceil(5*options_.psfSigmaPixels))+1;
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
    double radius = options_.radiusKm;
    if (radius == 0) {
        od::SpiceErrorModeGuard guard;
        SpiceDouble radii[3]; SpiceInt count;
        bodvrd_c(body.c_str(), "RADII", 3, &count, radii);
        od::throwIfSpiceFailed("Synthetic target PCK radii");
        require(count == 3, "Target needs three PCK radii or an explicit sphere radius.");
        radius = (radii[0]+radii[1]+radii[2])/3;
    }
    require(std::isfinite(radius) && radius > 0 && radius/range < 0.01,
            "Renderer requires a distant sphere (radius/range < 0.01).");
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
    truth.inertialToCamera = attitude.matrix();
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
    const Eigen::Matrix2d inverse = focal.inverse();
    const Eigen::Vector3d c = mid.centerCameraKm;
    // Local projected-area Jacobian. Full-phase Lambert disk integral is 2*pi*R^2/3.
    Eigen::Matrix<double,2,3> perspective;
    perspective << 1/c.z(),0,-c.x()/(c.z()*c.z()), 0,1/c.z(),-c.y()/(c.z()*c.z());
    const Eigen::Vector3d z = c.normalized(), x = z.unitOrthogonal(), y = z.cross(x);
    Eigen::Matrix<double,3,2> tangent; tangent.col(0)=x; tangent.col(1)=y;
    const double area = std::abs((focal*perspective*tangent).determinant());
    const double amplitude = options_.fullPhaseElectrons/(2*std::numbers::pi/3*mid.truth.radiusKm*mid.truth.radiusKm*area);
    require(std::isfinite(amplitude), "Invalid source radiance normalization.");
    const int sub = options_.subpixels, width = options_.size.width, height = options_.size.height;
    const int padding = int(std::ceil(5*options_.psfSigmaPixels))+1;
    // Render beyond the detector so off-frame light can still blur into it.
    image::Image high = image::Image::Zero((height+2*padding)*sub,(width+2*padding)*sub);
    const int samples = options_.exposureSeconds == 0 ? 1 : options_.exposureSamples;
    for (int sample = 0; sample < samples; ++sample) {
        const double offset = options_.exposureSeconds*((sample+0.5)/samples-0.5);
        const auto g = geometry(body,TdbEpoch(tdb.secondsPastJ2000+offset),observer,observerLabel);
        const double r2 = g.truth.radiusKm*g.truth.radiusKm;
        for (int l = 0; l < high.rows(); ++l) for (int s = 0; s < high.cols(); ++s) {
            const Eigen::Vector2d pixel{(s+0.5)/sub-padding-0.5-intrinsics.principalPoint.sample,
                                         (l+0.5)/sub-padding-0.5-intrinsics.principalPoint.line};
            const Eigen::Vector2d plane = inverse*pixel;
            const Eigen::Vector3d ray = Eigen::Vector3d(plane.x(),plane.y(),1).normalized();
            const double along = ray.dot(g.centerCameraKm);
            const double discriminant = r2-g.centerCameraKm.cross(ray).squaredNorm();
            if (along <= 0 || discriminant <= 0) continue;
            const double distance = along-std::sqrt(discriminant);
            const Eigen::Vector3d normal = (distance*ray-g.centerCameraKm)/g.truth.radiusKm;
            high(l,s) += amplitude*std::max(0.0,normal.dot(g.sunCamera))/(sub*sub*samples);
        }
    }
    const auto convolved = blur(high,options_.psfSigmaPixels*sub);
    image::Image source(height,width);
    for (int l = 0; l < height; ++l) for (int s = 0; s < width; ++s)
        source(l,s) = convolved.block((l+padding)*sub,(s+padding)*sub,sub,sub).sum();
    mid.truth.brightnessCenter = centerOfLight(source); mid.truth.seed = noise.seed;
    SyntheticOpticalFrame result;
    result.expectedSourceElectrons = source; result.truth = std::move(mid.truth);
    result.dn.resize(height,width); result.varianceDn2.resize(height,width);
    result.mask = image::PixelMask::Zero(height,width);
    std::mt19937_64 engine(noise.seed);
    std::normal_distribution<double> read(0,1);
    for (int l = 0; l < height; ++l) for (int s = 0; s < width; ++s) {
        const double expected = source(l,s)+noise.backgroundElectrons;
        require(std::isfinite(expected) && expected <= 1e12, "Unsupported detector electron count.");
        double electrons = expected;
        if (noise.shotNoise && expected > 0) electrons = double(std::poisson_distribution<std::uint64_t>(expected)(engine));
        if (noise.readNoiseElectrons > 0) electrons += noise.readNoiseElectrons*read(engine);
        double dn = noise.biasDn+electrons/noise.electronsPerDn;
        if (noise.quantize) dn = std::round(dn);
        result.mask(l,s) = dn >= noise.saturationDn || dn < 0;
        result.dn(l,s) = std::clamp(dn,0.0,noise.saturationDn);
        result.varianceDn2(l,s) = ((noise.shotNoise ? expected : 0)+noise.readNoiseElectrons*noise.readNoiseElectrons)/
            (noise.electronsPerDn*noise.electronsPerDn)+(noise.quantize ? 1.0/12 : 0);
    }
    return result;
}
} // namespace fd::opnav::synth
