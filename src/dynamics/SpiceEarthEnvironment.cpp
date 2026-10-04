#include "dynamics/SpiceEarthEnvironment.hpp"
#include "stations/StationCatalog.hpp"
#include "utils/CSPICE/SpiceError.hpp"
#include <SpiceUsr.h>
#include <cmath>

namespace fd::dynamics {
namespace {
double bodyConstant(const char* body, const char* item) {
    SpiceInt size;
    SpiceDouble values[3];
    bodvrd_c(body, item, 3, &size, values);
    od::throwIfSpiceFailed(std::string("Cannot read ") + body + " " + item);
    return values[0];
}
} // namespace

SpiceEarthEnvironment::SpiceEarthEnvironment(const std::filesystem::path& directory) {
    od::SpiceErrorModeGuard guard;
    try {
        for (const char* name : {"naif0012.tls", "pck00011.tpc", "de442.bsp", "gm_de440.tpc",
                                 "earthstns_itrf93_201023.bsp",
                                 "earth_1962_250826_2125_combined.bpc", "earth_assoc_itrf93.tf"}) {
            const auto path = std::filesystem::absolute(directory / name).string();
            furnsh_c(path.c_str());
            od::throwIfSpiceFailed("Cannot load " + path);
            kernels_.push_back(path);
        }
        earthMu_ = bodyConstant("EARTH", "GM");
        sunMu_ = bodyConstant("SUN", "GM");
        earthRadius_ = bodyConstant("EARTH", "RADII");
        sunRadius_ = bodyConstant("SUN", "RADII");
    } catch (...) {
        for (const auto& path : kernels_)
            unload_c(path.c_str());
        throw;
    }
}

SpiceEarthEnvironment::~SpiceEarthEnvironment() {
    od::SpiceErrorModeGuard guard;
    for (auto it = kernels_.rbegin(); it != kernels_.rend(); ++it)
        unload_c(it->c_str());
    reset_c();
}

double SpiceEarthEnvironment::epochTdb(const std::string& utc) const {
    od::SpiceErrorModeGuard guard;
    SpiceDouble epoch;
    str2et_c(utc.c_str(), &epoch);
    od::throwIfSpiceFailed("Cannot convert UTC start epoch");
    return epoch;
}

Eigen::Vector3d SpiceEarthEnvironment::sunPosition(double epoch) const {
    od::SpiceErrorModeGuard guard;
    SpiceDouble position[3], lightTime;
    spkpos_c("SUN", epoch, "J2000", "NONE", "EARTH", position, &lightTime);
    od::throwIfSpiceFailed("Cannot evaluate Earth-Sun geometry");
    return Eigen::Map<Eigen::Vector3d>(position);
}

CartesianState SpiceEarthEnvironment::stationState(const od::Station& station, double epoch) const {
    od::SpiceErrorModeGuard guard;
    SpiceDouble state[6], lightTime;
    const auto id = std::to_string(od::stationNaifIdFromName(station.name()));
    spkezr_c(id.c_str(), epoch, "J2000", "NONE", "EARTH", state, &lightTime);
    od::throwIfSpiceFailed("Cannot evaluate " + station.name() + " state");
    return {Eigen::Map<Eigen::Vector3d>(state), Eigen::Map<Eigen::Vector3d>(state + 3)};
}

Eigen::Vector3d SpiceEarthEnvironment::stationUp(const od::Station& station, double epoch) const {
    od::SpiceErrorModeGuard guard;
    SpiceDouble rotation[3][3];
    pxform_c("ITRF93", "J2000", epoch, rotation);
    od::throwIfSpiceFailed("Cannot transform geodetic station normal");
    const double latitude = station.latitudeRad(), longitude = station.longitudeRad();
    const Eigen::Vector3d up(std::cos(latitude) * std::cos(longitude),
                             std::cos(latitude) * std::sin(longitude), std::sin(latitude));
    return Eigen::Map<Eigen::Matrix<double, 3, 3, Eigen::RowMajor>>(&rotation[0][0]) * up;
}
} // namespace fd::dynamics
