#pragma once
#include "opnav/CameraModel.hpp"
#include "opnav/image/CentroidEstimator.hpp"
#include <optional>
#include <string>

namespace fd::opnav {

struct OpticalImage {
    std::string body;
    TdbEpoch time{0};
    double exposureSeconds{0};
    CameraAttitude attitude{Eigen::Matrix3d::Identity()}; // Recorded nominal rotation.
    Eigen::Matrix3d attitudeCovariance{Eigen::Matrix3d::Zero()}; // Camera-axis radians^2.
    Eigen::Vector2d smearPixels{Eigen::Vector2d::Zero()}; // Nominal end-minus-start displacement.
    std::optional<Eigen::Matrix3d> bodyToView; // Epoch-specific nominal orientation, for an ellipsoid.
    image::Image dn, varianceDn2;
    image::PixelMask mask;
};

struct CameraCalibration {
    CameraIntrinsics intrinsics;
    double psfSigmaPixels{0.8};
    double electronsPerDn{10}, readNoiseElectrons{3}, biasDn{100}, saturationDn{4095};
    bool digitized{true};
    Eigen::Matrix2d principalPointCovariance{Eigen::Matrix2d::Zero()};
    double relativeFocalSigma{0};
    double centerModelSigmaPixels{0}; // Explicit model-error allowance, not inferred accuracy.
};

struct BodyParameters {
    PixelCoordinates center;
    double radiusPixels{7}, phaseRadians{1.15}, sunAngleRadians{0.35};
    double amplitudeDn{1000}, backgroundDn{105};
    [[nodiscard]] Eigen::Matrix<double,7,1> vector() const;
    static BodyParameters fromVector(const Eigen::Matrix<double,7,1>&);
};

struct BodyModel {
    std::string body{"ENCELADUS"};
    // Sphere: ones. Ellipsoid: physical axis ratios with a declared body orientation.
    Eigen::Vector3d axes{Eigen::Vector3d::Ones()};
    // Orthogonal map to sample-right, line-down, observer-facing coordinates.
    // This view convention is left handed; determinant may be +1 or -1.
    Eigen::Matrix3d bodyToView{Eigen::Matrix3d::Identity()};
};

struct RenderedBody {
    image::Image dn;
    Eigen::Matrix<double,Eigen::Dynamic,7> jacobian; // Row-major pixel order, physical parameters.
};

class BodyRenderer {
public:
    explicit BodyRenderer(BodyModel model = {}, int quadratureOrder = 16);
    [[nodiscard]] RenderedBody render(ImageSize size, PixelCoordinates origin,
        const BodyParameters&, double psfSigma, Eigen::Vector2d smear = Eigen::Vector2d::Zero(),
        int exposureSamples = 1, bool derivatives = false) const;
    [[nodiscard]] const BodyModel& model() const noexcept { return model_; }
    [[nodiscard]] BodyRenderer reoriented(const Eigen::Matrix3d& orientation) const {
        auto model=model_;model.bodyToView=orientation;return BodyRenderer(model,order_);
    }
private:
    BodyModel model_;
    int order_;
    Eigen::VectorXd x_,y_,nx_,ny_,nz_,weights_;
};

struct BodyFitOptions {
    int windowRadius{18}, maxIterations{80}, exposureSamples{1};
    bool inferVarianceFromDn{true};
    double maximumReducedChiSquared{3};
};

struct BodyFitResult {
    image::FitStatus status{image::FitStatus::NonConverged};
    std::string body;
    TdbEpoch time{0};
    BodyParameters parameters;
    std::optional<image::PixelMeasurement> measurement; // Includes declared independent nuisance budgets.
    Eigen::Matrix2d formalCenterCovariance{Eigen::Matrix2d::Zero()};
    double reducedChiSquared{0};
    bool modelMismatch{false};
    int iterations{0};
    PixelCoordinates origin;
    image::Image modelDn, residualDn, varianceDn2;
    // Necessary fit diagnostics only; does not certify absolute measurement accuracy.
    [[nodiscard]] bool passesFitChecks() const {
        return status==image::FitStatus::Converged && measurement.has_value() && !modelMismatch;
    }
};

// No truth state, truth center, or SPICE pool access. Geometry/model registration
// is upstream; the body ID is checked and the fitted observable is timestamped.
class BodyEstimator {
public:
    BodyEstimator(BodyRenderer renderer, CameraCalibration calibration, BodyFitOptions options = {});
    [[nodiscard]] BodyFitResult fit(const OpticalImage&, std::optional<BodyParameters> initialGuess = std::nullopt) const;
private:
    BodyRenderer renderer_;
    CameraCalibration calibration_;
    BodyFitOptions options_;
};
} // namespace fd::opnav
