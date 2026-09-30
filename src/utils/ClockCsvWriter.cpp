#include "utils/ClockCsvWriter.hpp"

#include <fstream>
#include <iomanip>
#include <limits>
#include <locale>
#include <stdexcept>

namespace fd::utils {

void writeClockCsv(const std::filesystem::path& path, const clocks::ClockHistory& history,
    ClockCsvWriteMode mode) {
    if (mode != ClockCsvWriteMode::FailIfExists && mode != ClockCsvWriteMode::Overwrite)
        throw std::invalid_argument("Unknown clock CSV write mode.");
    if (mode == ClockCsvWriteMode::FailIfExists && std::filesystem::exists(path))
        throw std::runtime_error("Clock CSV already exists: " + path.string());
    std::ofstream stream(path, std::ios::out | std::ios::trunc);
    if (!stream) throw std::runtime_error("Cannot open clock CSV: " + path.string());
    stream.imbue(std::locale::classic());
    stream << std::setprecision(std::numeric_limits<double>::max_digits10);
    stream << "time_s,bias_s,fractional_frequency";
    if (history.includesCovariance()) stream << ",sigma_bias_s,sigma_fractional_frequency";
    stream << '\n';
    for (const auto& sample : history.samples()) {
        stream << sample.time_s << ',' << sample.bias_s << ',' << sample.fractional_frequency;
        if (history.includesCovariance())
            stream << ',' << sample.sigma_bias_s << ',' << sample.sigma_fractional_frequency;
        stream << '\n';
    }
    stream.flush();
    if (!stream) throw std::runtime_error("Failed writing clock CSV: " + path.string());
    stream.close();
    if (!stream) throw std::runtime_error("Failed closing clock CSV: " + path.string());
}

} // namespace fd::utils
