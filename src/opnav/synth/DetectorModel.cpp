#include "opnav/synth/DetectorModel.hpp"
#include <algorithm>
#include <cmath>
#include <random>
#include <stdexcept>
namespace fd::opnav::synth {
DetectorSample sampleDetector(const image::Image& source,OpticalDetectorNoise noise) {
    if(!source.allFinite() || (source.array()<0).any() || !std::isfinite(noise.electronsPerDn) || noise.electronsPerDn<1e-12 ||
       !std::isfinite(noise.backgroundElectrons) || noise.backgroundElectrons<0 || noise.backgroundElectrons>1e12 ||
       !std::isfinite(noise.readNoiseElectrons) || noise.readNoiseElectrons<0 || noise.readNoiseElectrons>1e6 ||
       !std::isfinite(noise.biasDn) || noise.biasDn<0 || !std::isfinite(noise.saturationDn) || noise.saturationDn<=noise.biasDn)
        throw std::invalid_argument("Invalid detector settings or expected source.");
    DetectorSample result;result.dn.resize(source.rows(),source.cols());result.varianceDn2.resize(source.rows(),source.cols());
    result.mask=image::PixelMask::Zero(source.rows(),source.cols());
    std::mt19937_64 engine(noise.seed);std::normal_distribution<double> read(0,1);
    for(int l=0;l<source.rows();++l)for(int s=0;s<source.cols();++s) {
        const double expected=source(l,s)+noise.backgroundElectrons;
        if(expected>1e12)throw std::invalid_argument("Unsupported detector electron count.");
        double electrons=expected;
        if(noise.shotNoise && expected>0)electrons=double(std::poisson_distribution<std::uint64_t>(expected)(engine));
        if(noise.readNoiseElectrons>0)electrons+=noise.readNoiseElectrons*read(engine);
        double dn=noise.biasDn+electrons/noise.electronsPerDn;if(noise.quantize)dn=std::round(dn);
        result.mask(l,s)=dn>=noise.saturationDn || dn<0;result.dn(l,s)=std::clamp(dn,0.,noise.saturationDn);
        result.varianceDn2(l,s)=((noise.shotNoise ? expected : 0)+noise.readNoiseElectrons*noise.readNoiseElectrons)/
            (noise.electronsPerDn*noise.electronsPerDn)+(noise.quantize ? 1./12 : 0);
    }
    return result;
}
} // namespace fd::opnav::synth
