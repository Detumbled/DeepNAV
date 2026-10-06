#include "dynamics/SpiceKernelSet.hpp"
#include "dynamics/SpiceStateProvider.hpp"
#include "observations/radiometric/PropagationCorrections.hpp"
#include "observations/radiometric/SpiceLightTime.hpp"
#include "simulation/CassiniRadioKernels.hpp"
#include "utils/CSPICE/SpiceError.hpp"
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <numbers>
#include <stdexcept>
#include <string>
#include <vector>

namespace radio = fd::observations::radiometric;
namespace {
struct Options {
    std::string utc{"2004-10-10T18:13:06.046"}, station{"auto"};
    int samples{3};
    double step{600}, count{60}, uplinkHz{7.175e9}, ramp{0}, delayUs{1}, zenithM{2.3};
    double maskDegrees{10};
    bool shapiro{true};
};
Options parse(int argc, char **argv) {
    Options o;
    for (int i = 1; i < argc; ++i) {
        const std::string key = argv[i];
        if (key == "--no-shapiro") {
            o.shapiro = false;
            continue;
        }
        if (key == "--help") {
            std::cout
                << "cassini_two_way_demo [--start-utc UTC] [--station auto|DSS-14|DSS-43|DSS-63]\n"
                   "  [--samples N] [--step-s S] [--count-s S] [--uplink-hz HZ]\n"
                   "  [--ramp-hz-s RATE] [--transponder-delay-us US] [--zenith-delay-m M]\n"
                   "  [--elevation-mask-deg DEG] [--no-shapiro]\n";
            std::exit(0);
        }
        if (++i >= argc)
            throw std::invalid_argument("Missing value for " + key);
        const std::string value = argv[i];
        if (key == "--start-utc")
            o.utc = value;
        else if (key == "--station")
            o.station = value;
        else {
            std::size_t consumed;
            const double x = std::stod(value, &consumed);
            if (consumed != value.size() || !std::isfinite(x))
                throw std::invalid_argument("Invalid numeric value for " + key);
            if (key == "--samples") {
                if (x < 1 || x > 10000 || x != std::floor(x))
                    throw std::invalid_argument("Samples must be an integer in [1,10000].");
                o.samples = static_cast<int>(x);
            } else if (key == "--step-s")
                o.step = x;
            else if (key == "--count-s")
                o.count = x;
            else if (key == "--uplink-hz")
                o.uplinkHz = x;
            else if (key == "--ramp-hz-s")
                o.ramp = x;
            else if (key == "--transponder-delay-us")
                o.delayUs = x;
            else if (key == "--zenith-delay-m")
                o.zenithM = x;
            else if (key == "--elevation-mask-deg")
                o.maskDegrees = x;
            else
                throw std::invalid_argument("Unknown option " + key);
        }
    }
    if (o.step <= 0 || o.count <= 0 || o.uplinkHz <= 0 || o.delayUs < 0 || o.zenithM < 0 ||
        o.maskDegrees <= 0 || o.maskDegrees >= 90)
        throw std::invalid_argument("Invalid interval, frequency, delay or elevation mask.");
    if (o.station != "auto" && o.station != "DSS-14" && o.station != "DSS-43" &&
        o.station != "DSS-63")
        throw std::invalid_argument("Unsupported DSN station.");
    return o;
}
std::string utc(double et) {
    char text[80];
    et2utc_c(et, "ISOC", 6, sizeof(text), text);
    od::throwIfSpiceFailed("Format radio event UTC");
    return text;
}
std::pair<double, double> elevations(const std::string &station, const radio::TwoWaySolution &s) {
    return {radio::spiceStationElevationRadians(station, s.uplink.emissionTdb,
                                                s.uplink.receiver.positionKm -
                                                    s.uplink.emitter.positionKm),
            radio::spiceStationElevationRadians(station, s.downlink.receptionTdb,
                                                s.downlink.emitter.positionKm -
                                                    s.downlink.receiver.positionKm)};
}
} // namespace
int main(int argc, char **argv) {
    try {
        const auto o = parse(argc, argv);
        od::SpiceErrorModeGuard errors("RETURN", "NONE");
        fd::dynamics::SpiceKernelSet kernels(fd::simulation::cassiniRadioKernels(
            std::filesystem::path(DEEPNAV_SOURCE_DIR) / "Kernels"));
        const radio::SpiceTdbClock clock;
        double start;
        str2et_c(o.utc.c_str(), &start);
        od::throwIfSpiceFailed("Parse reception UTC");
        fd::dynamics::SpiceStateProvider cassini("CASSINI"), sun("SUN");
        SpiceInt n;
        double sunGm;
        bodvrd_c("SUN", "GM", 1, &n, &sunGm);
        od::throwIfSpiceFailed("Solar GM");
        const radio::LinkConfig config{880.0 / 749.0, o.delayUs * 1e-6};
        const radio::TwoWayLink spice(radio::spiceReceptionLeg, config);
        const radio::TwoWayLink library(radio::libraryReceptionSolver(), config);
        const radio::FrequencyRamp ramp{start, o.uplinkHz, o.ramp};
        const double referenceHz = config.turnaroundRatio * o.uplinkHz;
        std::cout << std::fixed << std::setprecision(6)
                  << "Cassini coherent X/X two-way | ratio 880/749 = " << config.turnaroundRatio
                  << "\nUplink " << o.uplinkHz << " Hz; ramp " << o.ramp
                  << " Hz/s; constant receiver reference " << referenceHz << " Hz\n"
                  << "End-tagged count " << o.count << " TDB s; carrier clock TT (LSK TDB-TT).\n"
                  << "Illustrative non-dispersive transponder delay " << o.delayUs
                  << " us; troposphere zenith path " << o.zenithM << " m; solar Shapiro "
                  << (o.shapiro ? "ON" : "OFF") << "; elevation mask " << o.maskDegrees << " deg\n";
        double maxLtNs = 0, maxDopplerHz = 0, maxClockOffsetNs = 0;
        int valid = 0;
        for (int sample = 0; sample < o.samples; ++sample) {
            const double et = start + sample * o.step;
            std::string chosen;
            double best = -90;
            // Require visibility of both legs at both count endpoints. Auto-selection
            // changes only between counts; the station stays fixed during each count.
            for (const std::string name : {"DSS-14", "DSS-43", "DSS-63"}) {
                if (o.station != "auto" && o.station != name)
                    continue;
                fd::dynamics::SpiceStateProvider station(name);
                const auto a = elevations(name, spice.solve(et - o.count, station, cassini));
                const auto b = elevations(name, spice.solve(et, station, cassini));
                const double minimum =
                    std::min({a.first, a.second, b.first, b.second}) * 180 / std::numbers::pi;
                if (minimum > best) {
                    best = minimum;
                    chosen = name;
                }
            }
            if (chosen.empty() || best < o.maskDegrees) {
                std::cout << "\n"
                          << utc(et) << " SKIP: no station visible on both legs/count endpoints"
                          << " (best minimum elevation " << best << " deg)\n";
                continue;
            }
            fd::dynamics::SpiceStateProvider station(chosen);
            const auto a = spice.count(et, o.count, station, cassini, ramp, referenceHz, clock);
            const auto b = library.count(et, o.count, station, cassini, ramp, referenceHz, clock);
            for (double epoch : {et, et - o.count, a.end.uplink.emissionTdb}) {
                const double tt = unitim_c(epoch, "TDB", "TDT");
                od::throwIfSpiceFailed("Validate TDB-TT offset");
                const double error = std::abs((epoch - tt) - clock(epoch));
                maxClockOffsetNs = std::max(maxClockOffsetNs, error * 1e9);
                // unitim subtracts absolute epochs, so this reference is limited
                // by the double-precision epoch resolution (about 30 ns here).
                const double ulp = std::nextafter(epoch, epoch + 1) - epoch;
                if (error > 2 * ulp)
                    throw std::runtime_error("Smooth TDB-TT offset disagrees with CSPICE unitim.");
            }
            const auto delay = [&](const radio::LegSolution &leg, const auto &emitter,
                                   const auto &) {
                double seconds = o.shapiro ? radio::pointMassShapiroSeconds(leg, sun, sunGm) : 0;
                if (o.zenithM > 0) {
                    const bool uplink = &emitter == &station;
                    const auto direction = uplink
                                               ? leg.receiver.positionKm - leg.emitter.positionKm
                                               : leg.emitter.positionKm - leg.receiver.positionKm;
                    const double el = radio::spiceStationElevationRadians(
                        chosen, uplink ? leg.emissionTdb : leg.receptionTdb, direction);
                    if (el < o.maskDegrees * std::numbers::pi / 180)
                        throw std::runtime_error("Corrected link falls below station mask.");
                    // Simple dry, non-dispersive mapping; no measured meteorology.
                    seconds += radio::zenithPathDelaySeconds(o.zenithM, el);
                }
                return seconds;
            };
            const radio::TwoWayLink corrected(radio::libraryReceptionSolver({}, delay), config);
            const auto c = corrected.count(et, o.count, station, cassini, ramp, referenceHz, clock);
            const auto el = elevations(chosen, c.end);
            const double ltNs =
                std::max(
                    {std::abs(a.start.uplink.durationSeconds() - b.start.uplink.durationSeconds()),
                     std::abs(a.start.downlink.durationSeconds() -
                              b.start.downlink.durationSeconds()),
                     std::abs(a.end.uplink.durationSeconds() - b.end.uplink.durationSeconds()),
                     std::abs(a.end.downlink.durationSeconds() -
                              b.end.downlink.durationSeconds())}) *
                1e9;
            maxLtNs = std::max(maxLtNs, ltNs);
            maxDopplerHz = std::max(maxDopplerHz, std::abs(a.residualHz - b.residualHz));
            if (std::abs(a.residualHz - b.residualHz) > .005)
                throw std::runtime_error("CN/library Doppler benchmark exceeds 0.005 Hz.");
            if (ltNs > 10)
                throw std::runtime_error("CN/library benchmark exceeds 10 ns.");
            ++valid;
            std::cout << "\nReception " << utc(et) << " | " << chosen << " | elevation up/down "
                      << el.first * 180 / std::numbers::pi << "/"
                      << el.second * 180 / std::numbers::pi << " deg\n"
                      << "  Ground TX " << utc(c.end.uplink.emissionTdb) << "; Cassini RX "
                      << utc(c.end.spacecraftReceptionTdb) << "; Cassini TX "
                      << utc(c.end.spacecraftTransmissionTdb) << "\n"
                      << "  CN LT up/down " << a.end.uplink.durationSeconds() << "/"
                      << a.end.downlink.durationSeconds() << " s; RTLT " << a.end.roundTripSeconds()
                      << " s; c*RTLT/2 " << a.end.rangeKm() << " km\n"
                      << "  Library-CN max leg difference " << ltNs << " ns; range delta "
                      << (b.end.roundTripSeconds() - a.end.roundTripSeconds()) * 500 *
                             radio::speedOfLightKmPerSecond
                      << " m\n"
                      << "  Count Doppler received-reference: CN " << a.residualHz
                      << " Hz; library " << b.residualHz << " Hz; delta "
                      << b.residualHz - a.residualHz << " Hz\n"
                      << "  Corrected: RTLT " << c.end.roundTripSeconds() << " s; delay up/down "
                      << c.end.uplink.correctionSeconds * 1e6 << "/"
                      << c.end.downlink.correctionSeconds * 1e6 << " us; range correction "
                      << (c.end.roundTripSeconds() - b.end.roundTripSeconds()) * 500 *
                             radio::speedOfLightKmPerSecond
                      << " m\n"
                      << "  Corrected count Doppler " << c.residualHz << " Hz; correction "
                      << c.residualHz - b.residualHz << " Hz; mean carrier " << c.meanReceivedHz
                      << " Hz; count " << c.receivedCycles << " cycles\n";
        }
        std::cout
            << "\nBenchmark: " << valid << "/" << o.samples << " visible counts; max LT delta "
            << maxLtNs << " ns; max Doppler delta " << maxDopplerHz << " Hz.\n"
            << "TT offset check vs unitim: max " << maxClockOffsetNs
            << " ns (absolute-epoch reference resolution).\n"
            << "Synthetic observables, not recorded Cassini measurements. Range uses TDB RTLT.\n"
            << "No plasma/ionosphere, station hardware calibration, antenna offsets, ground clock\n"
            << "noise, spacecraft proper-time delay conversion or full relativistic Doppler.\n"
            << "Troposphere uses an illustrative zenith path; Shapiro is first-order solar only.\n";
        if (!valid)
            throw std::runtime_error("No visible counts; change reception time or station.");
        return 0;
    } catch (const std::exception &e) {
        std::cerr << "Cassini two-way demo: " << e.what() << '\n';
        return 1;
    }
}
