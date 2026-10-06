#include "dynamics/SpiceKernelSet.hpp"
#include "observations/radiometric/SpiceLightTime.hpp"
#include "utils/CSPICE/SpiceError.hpp"
#include <cmath>
#include <iostream>
namespace radio = fd::observations::radiometric;
int main() {
    try {
        od::SpiceErrorModeGuard guard("RETURN", "NONE");
        fd::dynamics::SpiceKernelSet kernels(
            {std::filesystem::path(DEEPNAV_SOURCE_DIR) / "Kernels/naif0012.tls"});
        const radio::SpiceTdbClock clock;
        for (double epoch : {0., 1.5e8, 8e8}) {
            const double tt = unitim_c(epoch, "TDB", "TDT");
            od::throwIfSpiceFailed("TT reference");
            const double ulp = std::nextafter(epoch, epoch + 1) - epoch;
            if (std::abs(epoch - tt - clock(epoch)) > std::max(2 * ulp, 1e-12))
                throw std::runtime_error("Cached clock disagrees with CSPICE TT");
        }
        const double epoch = 1.5e8, original = clock(epoch);
        double k;
        SpiceInt n;
        SpiceBoolean found;
        gdpool_c("DELTET/K", 0, 1, &n, &k, &found);
        od::throwIfSpiceFailed("Read original K");
        const double changed = 2 * k;
        pdpool_c("DELTET/K", 1, &changed);
        od::throwIfSpiceFailed("Change test clock kernel constant");
        const radio::SpiceTdbClock refreshed;
        if (clock(epoch) != original || std::abs(refreshed(epoch) - 2 * original) > 1e-15)
            throw std::runtime_error("Clock snapshot lifecycle failed");
        pdpool_c("DELTET/K", 1, &k);
        od::throwIfSpiceFailed("Restore clock constant");
        std::cout << "Cached LSK clock agrees with TT and reloads only on model recreation.\n";
        return 0;
    } catch (const std::exception &e) {
        std::cerr << e.what() << '\n';
        return 1;
    }
}
