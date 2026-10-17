#pragma once
#include "opnav/estimation/BodyEstimator.hpp"
#include <array>
#include <span>
#include <vector>

namespace fd::opnav {

enum class RejectionReason {
    FitFailed, InvalidInput, InvalidCovariance, PoorResiduals, InsufficientSignal,
    ExcessiveMasking, InsufficientBoundaryMargin, ExcessiveCenterUncertainty,
    AmbiguousSolution, InconsistentWithNavigation, MissingRequiredCheck
};
enum class RejectionCheck { FitValidity, Residuals, Signal, Masking, Boundary,
    Uncertainty, Ambiguity, Navigation, Count };
enum class CheckState { NotRun, Passed, Failed };

struct RejectionOptions {
    double maximumReducedChiSquared{3};
    double minimumPeakSignalToNoise{8}, minimumSourceSignalToNoise{25};
    double maximumMaskedFraction{.10};
    int minimumValidPixels{30};
    double psfSupportSigmas{3}, minimumBoundaryMarginPixels{1};
    double maximumFormalCenterSigmaPixels{.25}, maximumMeasurementCenterSigmaPixels{1};
    double ambiguityCenterSeparationPixels{.25}, ambiguityDeltaChiSquared{9};
    double maximumInnovationSquared{13.815510557964274}; // 99.9%, two dimensions.
    bool checkSignal{true}, checkMasking{true}, checkBoundary{true}, checkUncertainty{true};
    bool requireAmbiguityCheck{false}, requireNavigationCheck{false};
};

struct NavigationPrediction {
    std::string body;
    TdbEpoch time{0};
    PixelCoordinates center;
    Eigen::Matrix2d covariancePixels2{Eigen::Matrix2d::Zero()};
    // Cov(fitted measurement error, predicted measurement error). Innovation
    // covariance is R + P - C - C^T. Shared budgets must only be counted once.
    Eigen::Matrix2d measurementCrossCovariance{Eigen::Matrix2d::Zero()};
};

struct RejectionContext {
    const CameraCalibration* calibration{nullptr};
    const BodyRenderer* renderer{nullptr}; // Registered nominal shape; required for boundary/ambiguity checks.
    int exposureSamples{1};
    std::optional<NavigationPrediction> prediction;
    std::span<const BodyFitResult> alternatives;
};

struct RejectionDecision {
    bool evaluated{false};
    std::vector<RejectionReason> reasons;
    std::array<CheckState,static_cast<std::size_t>(RejectionCheck::Count)> checks{};
    std::optional<double> peakSignalToNoise, sourceSignalToNoise, maskedFraction;
    std::optional<double> boundaryMarginPixels, formalCenterSigmaPixels, measurementCenterSigmaPixels;
    std::optional<double> innovationSquared, ambiguousCenterSeparationPixels;
    int alternativesCompared{0};
    [[nodiscard]] bool accepted() const { return evaluated && reasons.empty(); }
};

[[nodiscard]] const char* rejectionReasonName(RejectionReason reason);

// Stateless, truth-free decision layer. Never refits, modifies a measurement,
// or certifies that accepted measurements are unbiased or mutually independent.
class EstimatorRejection {
public:
    explicit EstimatorRejection(RejectionOptions options = {});
    [[nodiscard]] RejectionDecision evaluate(const OpticalImage&, const BodyFitResult&,
        const RejectionContext& context = {}) const;
private:
    RejectionOptions options_;
};
} // namespace fd::opnav
