#pragma once
#include "observations/radiometric/TwoWayLink.hpp"

namespace fd::observations::radiometric {
struct VlbiSolution {
    LegSolution first, second;
    // Station-two reception minus station-one reception, coordinate TDB.
    // Both stations receive the same spacecraft emission event. Keep the
    // small difference directly rather than subtracting absolute epochs.
    double differentialSeconds{}, differentialRangeKm{}, residualSeconds{};
    std::size_t iterations{};
};
class VlbiLink {
  public:
    explicit VlbiLink(LegSolver firstLeg, DelayModel secondDelay = {},
                      double toleranceSeconds = 1e-12, std::size_t maxIterations = 16);
    [[nodiscard]] VlbiSolution solve(double firstReceptionTdb,
        const dynamics::StateProvider& spacecraft, const dynamics::StateProvider& stationOne,
        const dynamics::StateProvider& stationTwo) const;
  private:
    LegSolver firstLeg_;
    DelayModel secondDelay_;
    double tolerance_;
    std::size_t maxIterations_;
};
} // namespace fd::observations::radiometric
