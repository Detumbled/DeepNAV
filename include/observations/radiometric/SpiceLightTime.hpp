#pragma once
#include "observations/radiometric/TwoWayLink.hpp"
#include <string>

namespace fd::observations::radiometric {
// Requires SpiceStateProvider endpoints. CN, without stellar aberration or GR.
[[nodiscard]] LegSolution spiceReceptionLeg(double receptionTdb,
                                            const dynamics::StateProvider &emitter,
                                            const dynamics::StateProvider &receiver);
// Geodetic DSN elevation from the matching <station>_TOPO kernel frame.
[[nodiscard]] double spiceStationElevationRadians(const std::string &station, double epochTdb,
                                                  const Eigen::Vector3d &stationToTargetKm);
// Snapshot of LSK constants: create after loading kernels, recreate if the LSK
// changes. Evaluation never reads or mutates the process-global kernel pool.
class SpiceTdbClock {
  public:
    SpiceTdbClock();
    [[nodiscard]] double operator()(double epochTdb) const;

  private:
    double k_{}, eb_{}, meanAnomaly_{}, anomalyRate_{};
};
// Convenience one-shot query. For repeated counts, use a SpiceTdbClock instance.
// Smooth TDB-TT offset using the loaded LSK's DELTET periodic model. Unlike
// subtracting two absolute ET/TT epochs, preserves sub-microsecond differences.
[[nodiscard]] double spiceTdbMinusTt(double epochTdb);
} // namespace fd::observations::radiometric
