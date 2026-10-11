#include "observations/radiometric/VlbiLink.hpp"
#include <cmath>
#include <stdexcept>
#include <utility>

namespace fd::observations::radiometric {
namespace {
dynamics::CartesianState atOffset(const dynamics::StateProvider& provider, double reference,
                                 double offset) {
    const long double exact=static_cast<long double>(reference)+offset;
    const double rounded=static_cast<double>(exact);
    auto state=provider.stateAt(opnav::TdbEpoch{rounded});
    // Restore the sub-ULP station motion lost when adding a ~0.03 s delay to
    // a large epoch. Neglected acceleration over <one epoch ULP is negligible.
    state.positionKm+=state.velocityKmPerSec*(offset-(rounded-reference));
    if (!state.allFinite()) throw std::runtime_error("Nonfinite VLBI station state");
    return state;
}
}
VlbiLink::VlbiLink(LegSolver firstLeg, DelayModel secondDelay, double tolerance,
                   std::size_t maxIterations)
    :firstLeg_(std::move(firstLeg)),secondDelay_(std::move(secondDelay)),
     tolerance_(tolerance),maxIterations_(maxIterations) {
    if (!firstLeg_ || !std::isfinite(tolerance) || tolerance<=0 || maxIterations==0)
        throw std::invalid_argument("Invalid VLBI solver configuration");
}
VlbiSolution VlbiLink::solve(double reception, const dynamics::StateProvider& spacecraft,
    const dynamics::StateProvider& one, const dynamics::StateProvider& two) const {
    if (!std::isfinite(reception) || &one==&two || spacecraft.referenceFrame().empty() ||
        spacecraft.referenceOrigin().empty() || spacecraft.referenceFrame()!=one.referenceFrame() ||
        spacecraft.referenceFrame()!=two.referenceFrame() ||
        spacecraft.referenceOrigin()!=one.referenceOrigin() ||
        spacecraft.referenceOrigin()!=two.referenceOrigin())
        throw std::invalid_argument("VLBI requires distinct stations in a common inertial frame/origin");
    VlbiSolution result;result.first=firstLeg_(reception,spacecraft,one);
    if (!result.first.emitter.allFinite() || !result.first.receiver.allFinite() ||
        !std::isfinite(result.first.correctionSeconds)) throw std::runtime_error("Invalid first VLBI leg");
    result.second.emitter=result.first.emitter;
    result.second.emissionTdb=result.first.emissionTdb;
    using LongVector=Eigen::Matrix<long double,3,1>;
    const LongVector firstRay=result.first.emitter.positionKm.cast<long double>()-
                              result.first.receiver.positionKm.cast<long double>();
    const long double firstRange=firstRay.norm();
    if (!(firstRange>0)) throw std::runtime_error("Degenerate VLBI first leg");
    double offset=0;
    for (std::size_t iteration=1;iteration<=maxIterations_;++iteration) {
        result.second.receptionTdb=reception+offset;
        result.second.receiver=atOffset(two,reception,offset);
        const LongVector baseline=result.first.receiver.positionKm.cast<long double>()-
                                  result.second.receiver.positionKm.cast<long double>();
        const long double secondRange=(firstRay+baseline).norm();
        if (!(secondRange>0)) throw std::runtime_error("Degenerate VLBI second leg");
        result.second.geometricSeconds=static_cast<double>(secondRange/speedOfLightKmPerSecond);
        result.second.correctionSeconds=secondDelay_?secondDelay_(result.second,spacecraft,two):0;
        // Rationalized norm difference retains mm precision at interplanetary ranges.
        const long double difference=baseline.dot(2*firstRay+baseline)/(firstRange+secondRange);
        const double next=static_cast<double>(difference/speedOfLightKmPerSecond+
            (static_cast<long double>(result.second.correctionSeconds)-result.first.correctionSeconds));
        if (!std::isfinite(next)) throw std::runtime_error("Nonfinite VLBI delay");
        result.differentialSeconds=next;result.differentialRangeKm=next*speedOfLightKmPerSecond;
        result.residualSeconds=std::abs(next-offset);result.iterations=iteration;
        if (result.residualSeconds<=tolerance_) return result;
        offset=next;
    }
    throw std::runtime_error("VLBI common-wavefront iteration failed to converge");
}
} // namespace fd::observations::radiometric
