#pragma once
#include "opnav/estimation/BodyEstimator.hpp"
namespace fd::opnav {
// Upstream geometry adapter. Epoch is the nominal target emission epoch, based
// on navigation's light-time estimate. No spacecraft truth is queried here.
[[nodiscard]] BodyModel makeSpiceBodyModel(const std::string& body,TdbEpoch emissionEpoch,
    const CameraAttitude& nominalCamera,const CameraIntrinsics& nominalIntrinsics);
}
