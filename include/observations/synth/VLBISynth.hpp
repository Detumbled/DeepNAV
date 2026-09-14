#pragma once

#include "observations/synth/obs_synth.hpp"

#include <string>

namespace fd::observations::synth {

struct VLBIConfig {
    std::string target {"-31"};
    std::string stationOneName {"DSS-43"};
    std::string stationTwoName {"DSS-63"};
    std::string frame {"J2000"};
    std::string aberrationCorrection {"LT"};
    double minimumElevationRad {0.0};
};

// Represents a single VLBI observation sample, including the epoch, truth value, noise, observed value, and sigma.
struct VlbiObservation {
    std::string utc;
    std::string stationOneName;
    std::string stationTwoName;
    std::string stationOneNaif;
    std::string stationTwoNaif;
    double epochTdb {0.0};
    double delayTruthKm {0.0};
    double delayNoiseKm {0.0};
    double delayObservedKm {0.0};
    double delaySigmaKm {0.0};
};

class VLBISynth final : public SyntheticObservation {
public:
    explicit VLBISynth(VLBIConfig vlbi = {}, NoiseConfig noise = {});

    [[nodiscard]] std::vector<SyntheticObservationSample> generate(double startTdb,
                                                                   double endTdb,
                                                                   double stepSeconds) override;
    [[nodiscard]] std::vector<SyntheticObservationSample> generate(double startTdb,
                                                                   double endTdb,
                                                                   double stepSeconds,
                                                                   const TargetStateProvider& targetProvider);

    [[nodiscard]] std::vector<SyntheticObservationSample> generate(double startTdb,
                                                                   double endTdb,
                                                                   double stepSeconds,
                                                                   const std::string& stationOneName,
                                                                   const std::string& stationTwoName,
                                                                   const std::string& target);
    [[nodiscard]] std::vector<SyntheticObservationSample> generate(double startTdb,
                                                                   double endTdb,
                                                                   double stepSeconds,
                                                                   const std::string& stationOneName,
                                                                   const std::string& stationTwoName,
                                                                   const std::string& target,
                                                                   const TargetStateProvider& targetProvider);

    [[nodiscard]] const VLBIConfig& vlbiConfig() const noexcept;

private:
    [[nodiscard]] std::vector<SyntheticObservationSample> generateWithConfig(double startTdb,
                                                                             double endTdb,
                                                                             double stepSeconds,
                                                                             const VLBIConfig& vlbi);
    [[nodiscard]] std::vector<SyntheticObservationSample> generateWithConfig(
        double startTdb,
        double endTdb,
        double stepSeconds,
        const VLBIConfig& vlbi,
        const TargetStateProvider& targetProvider);

    VLBIConfig vlbi_;
};

} // namespace fd::observations::synth
