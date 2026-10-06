#include "dynamics/LinearStateProvider.hpp"
#include "observations/radiometric/TwoWayLink.hpp"
#include "observations/radiometric/PropagationCorrections.hpp"
#include <cmath>
#include <iostream>
#include <stdexcept>

namespace radio = fd::observations::radiometric;
namespace {
void near(double actual, double expected, double tolerance, const char* message) {
    if (!std::isfinite(actual) || std::abs(actual - expected) > tolerance)
        throw std::runtime_error(message);
}
fd::dynamics::LinearStateProvider linear(double range, double velocity = 0, double epoch = 0) {
    return {fd::opnav::TdbEpoch{epoch}, {{range, 0, 0}, {velocity, 0, 0}}, "J2000", "SSB"};
}
void stationaryAndRamp() {
    const auto station = linear(0), spacecraft = linear(3e6);
    const radio::LinkConfig config{1.25, 1e-6};
    const radio::TwoWayLink link(radio::libraryReceptionSolver(), config);
    const auto solution = link.solve(100, station, spacecraft);
    const double rtlt = 6e6 / radio::speedOfLightKmPerSecond + 1e-6;
    near(solution.roundTripSeconds(), rtlt, 1e-12, "Stationary round trip");
    near(solution.rangeKm(), 3e6 + radio::speedOfLightKmPerSecond * .5e-6,
         1e-8, "Transponder delay contribution");
    const double f = 7.175e9, reference = config.turnaroundRatio * f;
    auto count = link.count(100, 60, station, spacecraft, {0, f, 0}, reference);
    near(count.residualHz, 0, 1e-6, "Stationary coherent carrier");
    count = link.count(100, 60, station, spacecraft, {0, f, 2}, reference);
    near(count.residualHz, config.turnaroundRatio * 2 * (70 - rtlt), 1e-6,
         "Frequency ramp integrated at transmission epoch");
    // Constant clock offset cancels from intervals and ramp-relative epochs.
    const auto shifted = link.count(100, 60, station, spacecraft, {0, f, 2}, reference,
                                    [](double) { return 123.; });
    near(shifted.residualHz, count.residualHz, 1e-9, "Constant station clock bias cancellation");
}
void recedingAndLargeEpoch() {
    const double speed = 10, ratio = 2.75, f = 7.175e9;
    const radio::TwoWayLink link(radio::libraryReceptionSolver(), {ratio, 0});
    const double expected = ratio * f * (-2 * speed / (radio::speedOfLightKmPerSecond + speed));
    for (double epoch : {0., 1.5e8}) {
        const auto station = linear(0, 0, epoch), spacecraft = linear(3e6, speed, epoch);
        const auto count = link.count(epoch + 100, 60, station, spacecraft,
                                      {epoch, f, 0}, ratio * f);
        near(count.residualHz, expected, .001, "Analytic two-way receding Doppler");
        if (count.end.uplink.emissionTdb >= count.end.spacecraftReceptionTdb ||
            count.end.spacecraftTransmissionTdb >= count.end.downlink.receptionTdb)
            throw std::runtime_error("Radio event ordering");
    }
}
void correctedAndValidation() {
    const auto station = linear(0), spacecraft = linear(3e6), sun = linear(-1.5e8);
    const radio::TwoWayLink geometric(radio::libraryReceptionSolver(), {1, 0});
    const radio::TwoWayLink delayed(radio::libraryReceptionSolver({},
        [](const radio::LegSolution&, const auto&, const auto&) { return 1e-5; }), {1, 0});
    near(delayed.solve(100, station, spacecraft).roundTripSeconds() -
         geometric.solve(100, station, spacecraft).roundTripSeconds(), 2e-5, 1e-12,
         "Both legs receive correction");
    const auto leg = geometric.solve(100, station, spacecraft).downlink;
    const double gm = 1.32712440018e11;
    const double expected = 2 * gm / std::pow(radio::speedOfLightKmPerSecond, 3) *
                            std::log(1.53e8 / 1.5e8);
    near(radio::pointMassShapiroSeconds(leg, sun, gm), expected, 1e-15,
         "Radial point-mass Shapiro");
    near(radio::zenithPathDelaySeconds(2.3, std::asin(.5)),
         4.6 / (1000 * radio::speedOfLightKmPerSecond), 1e-20, "Zenith path mapping");
    bool rejected = false;
    try { (void)radio::TwoWayLink(radio::libraryReceptionSolver(), {1, -1}); }
    catch (const std::invalid_argument&) { rejected = true; }
    if (!rejected) throw std::runtime_error("Negative transponder delay accepted");
    const fd::dynamics::LinearStateProvider wrong(fd::opnav::TdbEpoch{0},
                                                 {{0,0,0},{0,0,0}}, "J2000", "EARTH");
    rejected = false;
    try { (void)geometric.solve(100, wrong, spacecraft); }
    catch (const std::invalid_argument&) { rejected = true; }
    if (!rejected) throw std::runtime_error("Mixed origins accepted");
}
} // namespace
int main() {
    try {
        stationaryAndRamp();
        recedingAndLargeEpoch();
        correctedAndValidation();
        std::cout << "Two-way event timing, coherent Doppler, ramp, clock cancellation, corrections "
                     "and validation passed.\n";
        return 0;
    } catch (const std::exception& e) { std::cerr << e.what() << '\n'; return 1; }
}
