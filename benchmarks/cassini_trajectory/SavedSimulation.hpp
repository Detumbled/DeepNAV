#pragma once
#include "CassiniForceModel.hpp"
#include "dynamics/EphemerisInterpolator.hpp"
#include "dynamics/StateProvider.hpp"
#include <filesystem>
#include <fstream>
#include <sstream>

namespace fd::simulation {
class SavedSimulation final : public fd::dynamics::StateProvider {
  od::EphemerisInterpolator history_;
  const fd::dynamics::StateProvider &origin_;

public:
  SavedSimulation(const std::filesystem::path &path,
                  const fd::simulation::CassiniForceModel &force,
                  const fd::dynamics::StateProvider &origin)
      : origin_(origin) {
    std::ifstream input(path);
    if (!input)
      throw std::runtime_error("Cannot read retained simulation CSV");
    std::string line;
    std::getline(input, line);
    if (line.find(
            "epoch_tdb,utc,model,x_km,y_km,z_km,vx_km_s,vy_km_s,vz_km_s") ==
        std::string::npos)
      throw std::runtime_error("Unexpected retained trajectory CSV schema");
    while (std::getline(input, line)) {
      std::stringstream row(line);
      std::vector<std::string> fields;
      std::string field;
      while (std::getline(row, field, ','))
        fields.push_back(field);
      if (fields.size() < 10)
        throw std::runtime_error("Malformed retained trajectory row");
      if (fields[3] != "multibody_j2_srp_thermal")
        continue;
      const double t = std::stod(fields[1]);
      fd::simulation::CassiniState x, dx;
      for (int j = 0; j < 6; ++j)
        x[j] = std::stod(fields[4 + j]);
      dx << x.tail<3>(),
          force.evaluate(t, x.head<3>(), x.tail<3>()).acceleration;
      history_.addNode(t, x, dx);
    }
    if (history_.history().size() < 2)
      throw std::runtime_error("Missing saved dynamical trajectory");
  }
  double startEpoch() const { return history_.history().front().tdb; }
  double endEpoch() const { return history_.history().back().tdb; }
  fd::simulation::CassiniState relative(double t) const {
    return history_.interpolate(t);
  }
  fd::dynamics::CartesianState stateAt(fd::opnav::TdbEpoch t) const override {
    const auto x = relative(t.secondsPastJ2000);
    const auto origin = origin_.stateAt(t);
    return {x.head<3>() + origin.positionKm,
            x.tail<3>() + origin.velocityKmPerSec};
  }
  std::string_view referenceFrame() const noexcept override { return "J2000"; }
  std::string_view referenceOrigin() const noexcept override { return "SSB"; }
};
} // namespace fd::simulation
