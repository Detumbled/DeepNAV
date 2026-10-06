#include "dynamics/LinearStateProvider.hpp"
#include "dynamics/SpiceKernelSet.hpp"
#include "dynamics/SpiceStateProvider.hpp"
#include "observations/radiometric/PropagationCorrections.hpp"
#include "observations/radiometric/SpiceLightTime.hpp"
#include "simulation/CassiniRadioKernels.hpp"
#include "utils/CSPICE/SpiceError.hpp"
#include <algorithm>
#include <chrono>
#include <iomanip>
#include <iostream>
#include <vector>
namespace radio = fd::observations::radiometric;
volatile double checksum;
template <class F> void measure(const char *name, int n, F operation) {
    for (int i = 0; i < 100; ++i)
        checksum = operation(i);
    std::vector<double> samples;
    for (int batch = 0; batch < 9; ++batch) {
        double sum = 0;
        const auto start = std::chrono::steady_clock::now();
        for (int i = 0; i < n; ++i)
            sum += operation(i);
        const auto end = std::chrono::steady_clock::now();
        checksum = sum;
        samples.push_back(std::chrono::duration<double, std::micro>(end - start).count() / n);
    }
    std::sort(samples.begin(), samples.end());
    std::cout << name << "," << std::setprecision(10) << samples[4] << "," << samples.front() << ","
              << samples.back() << '\n';
}
int main() {
    od::SpiceErrorModeGuard errors("RETURN", "NONE");
    constexpr double ratio = 880.0 / 749.0, f = 7.175e9;
    const fd::dynamics::LinearStateProvider station(fd::opnav::TdbEpoch{0}, {{0, 0, 0}, {0, 0, 0}},
                                                    "J2000", "SSB");
    const fd::dynamics::LinearStateProvider probe(fd::opnav::TdbEpoch{0}, {{3e6, 0, 0}, {10, 0, 0}},
                                                  "J2000", "SSB");
    const radio::TwoWayLink core(radio::libraryReceptionSolver(), {ratio, 1e-6});
    measure("numerical_count", 10000, [&](int i) {
        return core.count(100 + i * .1, 60, station, probe, {0, f, 0}, ratio * f).residualHz;
    });
    fd::dynamics::SpiceKernelSet kernels(
        fd::simulation::cassiniRadioKernels(std::filesystem::path(DEEPNAV_SOURCE_DIR) / "Kernels"));
    const radio::SpiceTdbClock clock;
    double epoch;
    str2et_c("2004-10-10T18:13:06.046", &epoch);
    od::throwIfSpiceFailed("epoch");
    fd::dynamics::SpiceStateProvider ground("DSS-14"), cassini("CASSINI"), sun("SUN");
    SpiceInt n;
    double gm;
    bodvrd_c("SUN", "GM", 1, &n, &gm);
    od::throwIfSpiceFailed("GM");
    const radio::TwoWayLink cn(radio::spiceReceptionLeg, {ratio, 1e-6});
    const radio::TwoWayLink library(radio::libraryReceptionSolver(), {ratio, 1e-6});
    const auto corrections = [&](const radio::LegSolution &leg, const auto &emitter, const auto &) {
        const bool up = &emitter == &ground;
        const auto direction = up ? leg.receiver.positionKm - leg.emitter.positionKm
                                  : leg.emitter.positionKm - leg.receiver.positionKm;
        const double elevation = radio::spiceStationElevationRadians(
            "DSS-14", up ? leg.emissionTdb : leg.receptionTdb, direction);
        return radio::pointMassShapiroSeconds(leg, sun, gm) +
               radio::zenithPathDelaySeconds(2.3, elevation);
    };
    const radio::TwoWayLink corrected(radio::libraryReceptionSolver({}, corrections),
                                      {ratio, 1e-6});
    for (const auto &pair :
         {std::pair{"spice_CN_count", &cn}, std::pair{"spice_library_count", &library},
          std::pair{"spice_corrected_count", &corrected}}) {
        measure(pair.first, 1000, [&](int i) {
            return pair.second
                ->count(epoch + i * .1, 60, ground, cassini, {epoch, f, 0}, ratio * f,
                        radio::spiceTdbMinusTt)
                .residualHz;
        });
    }
}
