#include "dynamics/SpiceKernelSet.hpp"
#include "utils/CSPICE/SpiceError.hpp"

namespace fd::dynamics {
SpiceKernelSet::SpiceKernelSet(const std::vector<std::filesystem::path> &paths) {
    od::SpiceErrorModeGuard guard;
    try {
        for (const auto &path : paths) {
            const auto absolute = std::filesystem::absolute(path).string();
            furnsh_c(absolute.c_str());
            od::throwIfSpiceFailed("Load " + absolute);
            loaded_.push_back(absolute);
        }
    } catch (...) {
        release();
        throw;
    }
}
SpiceKernelSet::~SpiceKernelSet() { release(); }
void SpiceKernelSet::release() noexcept {
    od::SpiceErrorModeGuard guard;
    for (auto it = loaded_.rbegin(); it != loaded_.rend(); ++it) {
        unload_c(it->c_str());
        if (failed_c())
            reset_c();
    }
}
} // namespace fd::dynamics
