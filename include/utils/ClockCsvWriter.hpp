#pragma once

#include "Clocks/ClockHistory.hpp"

#include <filesystem>

namespace fd::utils {

enum class ClockCsvWriteMode { FailIfExists, Overwrite };

// Export after simulation. Default refuses existing files; overwrite is explicit.
void writeClockCsv(const std::filesystem::path& path, const clocks::ClockHistory& history,
    ClockCsvWriteMode mode = ClockCsvWriteMode::FailIfExists);

} // namespace fd::utils
