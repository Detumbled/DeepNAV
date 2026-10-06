#pragma once
#include <filesystem>
#include <vector>

namespace fd::dynamics {
// Caller-selected kernels, loaded in precedence order. Never clears other kernels.
// Like CSPICE's kernel pool, construction/destruction must not run concurrently.
class SpiceKernelSet {
  public:
    explicit SpiceKernelSet(const std::vector<std::filesystem::path> &paths);
    ~SpiceKernelSet();
    SpiceKernelSet(const SpiceKernelSet &) = delete;
    SpiceKernelSet &operator=(const SpiceKernelSet &) = delete;

  private:
    void release() noexcept;
    std::vector<std::string> loaded_;
};
} // namespace fd::dynamics
