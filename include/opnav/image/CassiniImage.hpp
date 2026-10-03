#pragma once

#include "opnav/image/CircularGaussian.hpp"
#include <filesystem>
#include <string>

namespace fd::opnav::image {

struct CassiniImage {
    Image dn; // Original detector DN, without radiometric calibration.
    std::string imageMidTimeUtc;
    std::string targetDescription;
    double electronsPerDn;
    double exposureMilliseconds;
};

// Cassini PDS3 detached-label EDR: fixed records, single-band BSQ,
// big-endian signed 16-bit HALF samples. Checks PDS and VICAR layouts agree.
// Label defaults to the image path with extension .LBL.
[[nodiscard]] CassiniImage loadCassiniImage(
    const std::filesystem::path& imagePath,
    std::filesystem::path labelPath = {});

} // namespace fd::opnav::image
