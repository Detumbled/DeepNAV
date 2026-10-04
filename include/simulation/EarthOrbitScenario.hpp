#pragma once

#include "Clocks/Types.hpp"
#include "perturbations/ForceModel.hpp"
#include <string>

namespace fd::dynamics {
class SpiceEarthEnvironment;
}

namespace fd::simulation {

inline constexpr double kNominalEarthJ2 = 1.08262668e-3;

struct EarthOrbitModel {
    double srpCr{1.3}, srpAreaM2{20}, srpMassKg{1000};
    fd::clocks::ClockParameters clock;
};

// Changes estimator assumptions relative to fixed truth; scales multiply parameters, not draws.
struct ModelMismatch {
    double srpCrScale{1}, srpAreaMassScale{1};
    double clockDriftOffsetPerDay{0}, clockNoiseScale{1};
};

[[nodiscard]] EarthOrbitModel nominalEarthOrbitModel();
[[nodiscard]] ModelMismatch namedModelMismatch(const std::string& name);
[[nodiscard]] EarthOrbitModel estimatorModel(const EarthOrbitModel& truth,
                                             const ModelMismatch& mismatch);
[[nodiscard]] fd::perturbations::AccelerationFunction
earthOrbitForces(const fd::dynamics::SpiceEarthEnvironment& environment, double startEpoch,
                 const EarthOrbitModel& model, bool twoBodyOnly = false);

} // namespace fd::simulation
