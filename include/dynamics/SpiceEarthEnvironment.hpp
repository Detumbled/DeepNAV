#pragma once

#include "dynamics/CartesianState.hpp"
#include "stations/Stations.hpp"
#include <Eigen/Core>
#include <filesystem>
#include <string>
#include <vector>

namespace fd::dynamics {

// Loads the demo's local kernels; CSPICE's kernel pool is process-global, not
// thread-safe.
class SpiceEarthEnvironment {
  public:
    explicit SpiceEarthEnvironment(const std::filesystem::path& kernelDirectory);
    ~SpiceEarthEnvironment();
    SpiceEarthEnvironment(const SpiceEarthEnvironment&) = delete;
    SpiceEarthEnvironment& operator=(const SpiceEarthEnvironment&) = delete;
    [[nodiscard]] double epochTdb(const std::string& utc) const;
    [[nodiscard]] Eigen::Vector3d sunPosition(double epochTdb) const;
    [[nodiscard]] CartesianState stationState(const od::Station& station, double epochTdb) const;
    [[nodiscard]] Eigen::Vector3d stationUp(const od::Station& station, double epochTdb) const;
    [[nodiscard]] double earthMu() const noexcept { return earthMu_; }
    [[nodiscard]] double sunMu() const noexcept { return sunMu_; }
    [[nodiscard]] double earthRadius() const noexcept { return earthRadius_; }
    [[nodiscard]] double sunRadius() const noexcept { return sunRadius_; }

  private:
    std::vector<std::string> kernels_;
    double earthMu_{}, sunMu_{}, earthRadius_{}, sunRadius_{};
};

} // namespace fd::dynamics
