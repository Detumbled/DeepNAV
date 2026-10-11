#include "observations/synth/VLBISynth.hpp"

#include "stations/StationCatalog.hpp"
#include "perturbations/Shapiro.hpp"
#include "observations/radiometric/VlbiLink.hpp"
#include "observations/radiometric/PropagationCorrections.hpp"
#include "utils/CSPICE/SpiceError.hpp"
#include <SpiceUsr.h>

#include <cmath>
#include <stdexcept>
#include <utility>

namespace fd::observations::synth {
namespace {

[[nodiscard]] GeometryConfig makePrimaryGeometry(const VLBIConfig& vlbi) {
    GeometryConfig geometry;
    geometry.target = vlbi.target;
    geometry.stationName = vlbi.stationOneName;
    geometry.frame = vlbi.frame;
    geometry.aberrationCorrection = vlbi.aberrationCorrection;
    return geometry;
}

void validateVLBIConfig(const VLBIConfig& vlbi) {
    if (vlbi.target.empty()) {
        throw std::invalid_argument("VLBI synthetic observation target cannot be empty.");
    }
    if (vlbi.stationOneName.empty() || vlbi.stationTwoName.empty()) {
        throw std::invalid_argument("VLBI synthetic observation stations cannot be empty.");
    }
    if (vlbi.stationOneName == vlbi.stationTwoName) {
        throw std::invalid_argument("VLBI synthetic observation requires two distinct stations.");
    }
    if (vlbi.frame.empty()) {
        throw std::invalid_argument("VLBI synthetic observation frame cannot be empty.");
    }
    if (vlbi.aberrationCorrection.empty()) {
        throw std::invalid_argument("VLBI synthetic observation aberration correction cannot be empty.");
    }
    if (!std::isfinite(vlbi.minimumElevationRad)) {
        throw std::invalid_argument("VLBI synthetic observation elevation mask must be finite.");
    }
}

[[nodiscard]] od::Station stationFromKernel(const std::string& stationName, double epochTdb) {
    const int station_naif_id = od::stationNaifIdFromName(stationName);
    return od::buildStationFromKernel(stationName, station_naif_id, epochTdb);
}

// Adapt the legacy Sun-relative provider to the shared inertial SSB link.
class Endpoint final:public fd::dynamics::StateProvider {
    std::string target_,frame_;const TargetStateProvider* injected_;
public:
    Endpoint(std::string target,std::string frame,const TargetStateProvider* injected=nullptr):target_(std::move(target)),frame_(std::move(frame)),injected_(injected) {}
    fd::dynamics::CartesianState stateAt(fd::opnav::TdbEpoch epoch) const override {
        double state[6],lt;spkezr_c((injected_?"SUN":target_.c_str()),epoch.secondsPastJ2000,frame_.c_str(),"NONE","SSB",state,&lt);od::throwIfSpiceFailed("VLBI endpoint state");
        TargetState6 x;for (int j=0;j<6;++j) x[j]=state[j];if (injected_) x+=injected_->stateAt(epoch.secondsPastJ2000);
        return {x.head<3>(),x.tail<3>()};
    }
    std::string_view referenceFrame() const noexcept override {return frame_;}
    std::string_view referenceOrigin() const noexcept override {return "SSB";}
};

} // namespace

VLBISynth::VLBISynth(VLBIConfig vlbi, NoiseConfig noise)
    : SyntheticObservation(makePrimaryGeometry(vlbi), noise),
      vlbi_(std::move(vlbi)) {
    validateVLBIConfig(vlbi_);
}

std::vector<SyntheticObservationSample> VLBISynth::generate(double startTdb,
                                                            double endTdb,
                                                            double stepSeconds) {
    return generateWithConfig(startTdb, endTdb, stepSeconds, vlbi_);
}

std::vector<SyntheticObservationSample> VLBISynth::generate(
    double startTdb,
    double endTdb,
    double stepSeconds,
    const TargetStateProvider& targetProvider) {
    return generateWithConfig(startTdb, endTdb, stepSeconds, vlbi_, targetProvider);
}

std::vector<SyntheticObservationSample> VLBISynth::generate(double startTdb,
                                                            double endTdb,
                                                            double stepSeconds,
                                                            const std::string& stationOneName,
                                                            const std::string& stationTwoName,
                                                            const std::string& target) {
    VLBIConfig vlbi = vlbi_;
    vlbi.stationOneName = stationOneName;
    vlbi.stationTwoName = stationTwoName;
    vlbi.target = target;
    validateVLBIConfig(vlbi);
    return generateWithConfig(startTdb, endTdb, stepSeconds, vlbi);
}

std::vector<SyntheticObservationSample> VLBISynth::generate(
    double startTdb,
    double endTdb,
    double stepSeconds,
    const std::string& stationOneName,
    const std::string& stationTwoName,
    const std::string& target,
    const TargetStateProvider& targetProvider) {
    VLBIConfig vlbi = vlbi_;
    vlbi.stationOneName = stationOneName;
    vlbi.stationTwoName = stationTwoName;
    vlbi.target = target;
    validateVLBIConfig(vlbi);
    return generateWithConfig(startTdb, endTdb, stepSeconds, vlbi, targetProvider);
}

const VLBIConfig& VLBISynth::vlbiConfig() const noexcept {
    return vlbi_;
}

std::vector<SyntheticObservationSample> VLBISynth::generateWithConfig(double startTdb,
                                                                      double endTdb,
                                                                      double stepSeconds,
                                                                      const VLBIConfig& vlbi) {
    const SpiceTargetStateProvider target_provider(vlbi.target, vlbi.frame, "SUN");
    return generateWithConfig(startTdb, endTdb, stepSeconds, vlbi, target_provider);
}

std::vector<SyntheticObservationSample> VLBISynth::generateWithConfig(
    double startTdb,
    double endTdb,
    double stepSeconds,
    const VLBIConfig& vlbi,
    const TargetStateProvider& targetProvider) {
    const std::vector<double> epochs = makeEpochGrid(startTdb, endTdb, stepSeconds);
    const od::Station station_one = stationFromKernel(vlbi.stationOneName, startTdb);
    const od::Station station_two = stationFromKernel(vlbi.stationTwoName, startTdb);

    std::vector<SyntheticObservationSample> samples;
    samples.reserve(epochs.size());

    const Endpoint probe(vlbi.target,vlbi.frame,&targetProvider),one(vlbi.stationOneName,vlbi.frame),two(vlbi.stationTwoName,vlbi.frame),sun("SUN",vlbi.frame);
    namespace radio=fd::observations::radiometric;
    const radio::DelayModel delay=[&](const auto& leg,const auto&,const auto&) {return radio::pointMassShapiroSeconds(leg,sun,fd::perturbations::kSunGravitationalParameterKm3PerSec2);};
    const radio::VlbiLink link(radio::libraryReceptionSolver({},delay),delay);
    for (const double epoch : epochs) {
        const bool station_one_visible =
            targetElevationRadFor(vlbi.stationOneName, station_one, epoch, targetProvider)
            >= vlbi.minimumElevationRad;
        const bool station_two_visible =
            targetElevationRadFor(vlbi.stationTwoName, station_two, epoch, targetProvider)
            >= vlbi.minimumElevationRad;
        if (!station_one_visible || !station_two_visible) {
            continue;
        }

        const double truth = link.solve(epoch,probe,one,two).differentialRangeKm;
        const double noise = drawNoise();

        samples.push_back(SyntheticObservationSample{
            epoch,
            truth,
            noise,
            truth + noise,
            noiseConfig().sigma
        });
    }

    return samples;
}

} // namespace fd::observations::synth
