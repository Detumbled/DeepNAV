#pragma once
#include <filesystem>
#include <vector>
namespace fd::simulation {
// Mission data choices belong to the demo, not the generic radio interfaces.
inline std::vector<std::filesystem::path> cassiniRadioKernels(const std::filesystem::path &root) {
    std::vector<std::filesystem::path> paths;
    for (const char *name : {"naif0012.tls", "pck00011.tpc", "gm_de440.tpc", "de442.bsp",
                             "earthstns_itrf93_201023.bsp", "earth_1962_250826_2125_combined.bpc",
                             "earth_assoc_itrf93.tf", "Cassini/fk/earth_topo_201023.tf",
                             "Cassini/spk/200128RU_SCPSE_04251_04303.bsp"})
        paths.push_back(root / name);
    return paths;
}
} // namespace fd::simulation
