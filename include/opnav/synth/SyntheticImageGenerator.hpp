#pragma once

#include "opnav/CameraModel.hpp"
#include "opnav/image/CentroidEstimator.hpp"
#include "dynamics/StateProvider.hpp"
#include "opnav/synth/DetectorModel.hpp"
#include <functional>
#include <optional>
#include <string>

namespace fd::opnav::synth {

struct SyntheticImageOptions {
    ImageSize size{128, 128};
    double radiusKm{0}; // Zero: use mean of target PCK radii. Otherwise declared sphere radius.
    double fullPhaseElectrons{1e6}; // Lambert-sphere-equivalent budget; alternate surfaces are not renormalized.
    double psfSigmaPixels{1};
    double exposureSeconds{0.06};
    std::optional<double> psfSigmaLinePixels; // Elliptical truth PSF, sample sigma above.
    bool usePckEllipsoid{false};
    std::optional<Eigen::Vector3d> bodyAxesKm;
    double albedoContrast{0}; // Body-fixed smooth pattern, amplitude in [0,1).
    double lommelSeeligerMix{0}; // Blend with Lambert, in [0,1].
    double referenceRangeKm{0},referenceExposureSeconds{0.06}; // Optional physical count scaling.
    int exposureSamples{5}; // Uniform shutter, midpoint quadrature.
    int subpixels{4}; // Rendering and PSF convolution before detector pixel integration.
    // Optional synthetic lighting for controlled phase tests. Otherwise use SPK Sun.
    std::optional<Eigen::Vector3d> sunDirectionCamera;
};

struct OpticalPixelNoise {
    double sigmaPixels{0.1}; // Independent sample/line Gaussian noise.
    std::uint64_t seed{1};
};

struct OpticalFrameTruth {
    std::string target, observer;
    TdbEpoch epoch{0}; // Exposure midpoint, TDB seconds past J2000.
    PixelCoordinates geometricCenter;
    std::optional<PixelCoordinates> brightnessCenter; // Noiseless source after PSF and detector crop.
    double rangeKm{0}, lightTimeSeconds{0}, phaseAngleRadians{0}, radiusKm{0};
    Eigen::Vector3d observerPositionKm, observerVelocityKmPerSecond;
    Eigen::Vector3d bodyAxesKm, sunDirectionCamera;
    Eigen::Matrix3d bodyToCamera;
    Eigen::Matrix3d inertialToCamera;
    std::uint64_t seed{0};
};

struct SyntheticOpticalFrame {
    image::Image dn, varianceDn2;
    image::PixelMask mask; // Pixels clipped at either detector bound.
    image::Image expectedSourceElectrons; // Before noise/bias/background; includes exposure and PSF.
    OpticalFrameTruth truth; // Evaluation only; never pass to centroid estimation as a measurement.
};

struct SyntheticOpticalMeasurement {
    image::PixelMeasurement measurement;
    OpticalFrameTruth truth;
};

// Caller loads SPICE kernels and supplies truth attitude. No kernel-pool mutation.
// One distant body: sphere or triaxial ellipsoid, Lambert/LS blend and optional
// body-fixed albedo. Supports inverse camera distortion and an axis-aligned
// elliptical PSF. No terrain, stars, or mutually occulting bodies.
class SyntheticImageGenerator {
public:
    using AttitudeProvider = std::function<CameraAttitude(TdbEpoch)>;
    SyntheticImageGenerator(CameraModel camera, AttitudeProvider attitude,
                            SyntheticImageOptions options = {});
    [[nodiscard]] SyntheticOpticalFrame createImage(const std::string& body, TdbEpoch tdb,
        OpticalDetectorNoise noise, const std::string& observer) const;
    [[nodiscard]] SyntheticOpticalFrame create(const std::string& body, TdbEpoch tdb,
        OpticalDetectorNoise noise, const std::string& observer) const {
        return createImage(body, tdb, noise, observer);
    }
    // Geometric-center angular observation: does not simulate brightness-center bias.
    [[nodiscard]] SyntheticOpticalMeasurement createMeasurement(const std::string& body,
        TdbEpoch tdb, OpticalPixelNoise noise, const std::string& observer) const;
    // Inject simulated truth in geometric J2000/SSB coordinates. The provider is
    // queried across the exposure and is never replaced by a spacecraft SPK query.
    [[nodiscard]] SyntheticOpticalFrame createImage(const std::string& body, TdbEpoch tdb,
        OpticalDetectorNoise noise, const dynamics::StateProvider& observer,
        const std::string& observerLabel = "SIMULATED") const;
    [[nodiscard]] SyntheticOpticalFrame create(const std::string& body, TdbEpoch tdb,
        OpticalDetectorNoise noise, const dynamics::StateProvider& observer) const {
        return createImage(body,tdb,noise,observer);
    }
    [[nodiscard]] SyntheticOpticalMeasurement createMeasurement(const std::string& body,
        TdbEpoch tdb, OpticalPixelNoise noise, const dynamics::StateProvider& observer,
        const std::string& observerLabel = "SIMULATED") const;
private:
    struct Geometry;
    [[nodiscard]] Geometry geometry(const std::string& body, TdbEpoch epoch,
        const dynamics::StateProvider& observer, const std::string& observerLabel) const;
    CameraModel camera_;
    AttitudeProvider attitude_;
    SyntheticImageOptions options_;
};
} // namespace fd::opnav::synth
