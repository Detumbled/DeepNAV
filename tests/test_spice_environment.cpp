#include "dynamics/SpiceEarthEnvironment.hpp"
#include "stations/StationCatalog.hpp"
#include "utils/CSPICE/SpiceError.hpp"
#include <SpiceUsr.h>
#include <cmath>
#include <iostream>
#include <stdexcept>

namespace {
void require(bool condition, const char* message) {
    if (!condition)
        throw std::runtime_error(message);
}
} // namespace

int main() {
    try {
        const fd::dynamics::SpiceEarthEnvironment environment(
            std::filesystem::path(DEEPNAV_SOURCE_DIR) / "Kernels");
        const double epoch = environment.epochTdb("2024-03-20T00:00:00");
        require(environment.earthMu() > 398000 && environment.sunRadius() > 690000,
                "Kernel constants are invalid.");
        require((environment.sunPosition(epoch + 28800) - environment.sunPosition(epoch)).norm() >
                    100000,
                "Solar ephemeris is not advancing.");
        for (const char* name : {"DSS-43", "DSS-63", "DSS-14"}) {
            const auto station =
                od::buildStationFromKernel(name, od::stationNaifIdFromName(name), epoch);
            const auto state = environment.stationState(station, epoch);
            require(state.positionKm.norm() > 6300 && state.positionKm.norm() < 6400,
                    "Station position has invalid units/origin.");
            const Eigen::Vector3d derivative =
                (environment.stationState(station, epoch + 1).positionKm -
                 environment.stationState(station, epoch - 1).positionKm) /
                2;
            require((derivative - state.velocityKmPerSec).norm() < 1e-8,
                    "Station velocity disagrees with position derivative.");
            const auto up = environment.stationUp(station, epoch);
            require(std::abs(up.norm() - 1) < 1e-12 && up.dot(state.positionKm.normalized()) > .999,
                    "Geodetic station normal/frame is invalid.");
            // Independent ITRF93 -> J2000 full state transform checks frame/origin
            // consistency.
            od::SpiceErrorModeGuard guard;
            SpiceDouble fixed[6], transformed[6], transform[6][6], lightTime;
            const auto id = std::to_string(od::stationNaifIdFromName(name));
            spkezr_c(id.c_str(), epoch, "ITRF93", "NONE", "EARTH", fixed, &lightTime);
            sxform_c("ITRF93", "J2000", epoch, transform);
            mxvg_c(transform, fixed, 6, 6, transformed);
            od::throwIfSpiceFailed("Station state transform test");
            require(
                (state.positionKm - Eigen::Map<Eigen::Vector3d>(transformed)).norm() < 1e-8 &&
                    (state.velocityKmPerSec - Eigen::Map<Eigen::Vector3d>(transformed + 3)).norm() <
                        1e-12,
                "Station state disagrees with Earth-fixed frame transform.");
        }
        bool rejected = false;
        try {
            od::SpiceErrorModeGuard guard("RETURN", "NONE");
            (void)environment.stationState(od::buildStationFromKernel("DSS-14", 399014, epoch),
                                           environment.epochTdb("1800-01-01 UTC"));
        } catch (const std::runtime_error&) {
            rejected = true;
        }
        require(rejected, "Missing kernel coverage must fail explicitly.");
        std::cout << "SPICE station states, velocity, normals, epoch and coverage "
                     "checks passed.\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "FAIL: " << error.what() << '\n';
        return 1;
    }
}
