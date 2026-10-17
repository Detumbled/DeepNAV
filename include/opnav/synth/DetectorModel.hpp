#pragma once
#include "opnav/image/CircularGaussian.hpp"

namespace fd::opnav::synth {
struct OpticalDetectorNoise {
    double electronsPerDn{1},backgroundElectrons{20},readNoiseElectrons{3},biasDn{0},saturationDn{4095};
    bool shotNoise{true},quantize{true};
    std::uint64_t seed{1};
};
struct DetectorSample {image::Image dn,varianceDn2;image::PixelMask mask;};
// Independent sampling permits Monte Carlo reuse of the same noiseless source.
[[nodiscard]] DetectorSample sampleDetector(const image::Image& sourceElectrons,OpticalDetectorNoise noise);
} // namespace fd::opnav::synth
