#include "opnav/image/CassiniImage.hpp"
#include <cmath>
#include <fstream>
#include <regex>
#include <stdexcept>
#include <vector>

namespace fd::opnav::image {
namespace {
void check(bool condition, const char* message) {
    if (!condition) throw std::runtime_error(message);
}
std::string field(const std::string& text, const std::string& key) {
    std::smatch match;
    const std::regex pattern("(^|[\\r\\n])\\s*" + key + "\\s*=\\s*([^\\r\\n]+)");
    check(std::regex_search(text, match, pattern), ("Missing PDS field: " + key).c_str());
    std::string value = match[2];
    const auto end = value.find_last_not_of(" \t");
    value.resize(end + 1);
    if (value.size() >= 2 && value.front() == '"' && value.back() == '"')
        value = value.substr(1, value.size()-2);
    return value;
}
std::size_t number(const std::string& value) {
    check(!value.empty() && value.find_first_not_of("0123456789") == std::string::npos,
          "Invalid nonnegative integer in image layout.");
    return std::stoull(value);
}
std::string vicar(const std::string& header, const std::string& key) {
    std::smatch match;
    check(std::regex_search(header, match,
          std::regex("(^|\\s)" + key + "=\\s*('([^']*)'|([0-9]+))")),
          ("Missing VICAR field: " + key).c_str());
    return match[3].matched ? match[3].str() : match[4].str();
}
}

CassiniImage loadCassiniImage(const std::filesystem::path& imagePath,
                             std::filesystem::path labelPath) {
    if (labelPath.empty()) { labelPath = imagePath; labelPath.replace_extension(".LBL"); }
    std::ifstream labelFile(labelPath);
    check(bool(labelFile), "Cannot open Cassini .LBL file.");
    const std::string label((std::istreambuf_iterator<char>(labelFile)), {});
    check(field(label, "PDS_VERSION_ID") == "PDS3" &&
          field(label, "RECORD_TYPE") == "FIXED_LENGTH", "Unsupported PDS image format.");
    std::smatch object;
    check(std::regex_search(label, object,
          std::regex("(^|[\\r\\n])\\s*OBJECT\\s*=\\s*IMAGE\\s*[\\r\\n]([\\s\\S]*?)END_OBJECT\\s*=\\s*IMAGE")),
          "Missing PDS IMAGE object.");
    const std::string layout = object[2];
    check(field(layout, "SAMPLE_TYPE") == "SUN_INTEGER" &&
          number(field(layout, "SAMPLE_BITS")) == 16, "Expected big-endian signed 16-bit pixels.");
    const auto lines = number(field(layout, "LINES"));
    const auto samples = number(field(layout, "LINE_SAMPLES"));
    const auto prefix = number(field(layout, "LINE_PREFIX_BYTES"));
    const auto record = number(field(label, "RECORD_BYTES"));
    const auto records = number(field(label, "FILE_RECORDS"));
    const auto fileSize = std::filesystem::file_size(imagePath);
    check(lines > 0 && samples > 0 && record > 0 && records > 0 &&
          samples <= record/2 && prefix <= record && samples*2 == record-prefix,
          "Invalid PDS record layout.");
    check(fileSize % record == 0 && fileSize/record == records, "IMG size does not match PDS label.");
    std::smatch pointer;
    const auto imagePointer = field(label, "\\^IMAGE");
    check(std::regex_match(imagePointer, pointer,
          std::regex("\\(\\s*\"([^\"]+)\"\\s*,\\s*([0-9]+)\\s*\\)")),
          "Expected detached record-based PDS IMAGE pointer.");
    check(pointer[1].str() == imagePath.filename().string(), "PDS IMAGE filename does not match IMG.");
    const auto firstRecord = number(pointer[2]);
    check(firstRecord > 0 && firstRecord <= records && lines == records-firstRecord+1,
          "Invalid PDS IMAGE offset or row count.");
    check(number(field(label, "MISSING_LINES")) == 0 &&
          field(label, "MISSING_PACKET_FLAG") == "NO", "Incomplete Cassini image is unsupported.");
    check(field(label, "DATA_CONVERSION_TYPE") == "12BIT", "Only direct 12BIT DN conversion is supported.");

    std::ifstream input(imagePath, std::ios::binary);
    check(bool(input), "Cannot open Cassini IMG file.");
    char lead[64]{};
    input.read(lead, sizeof lead);
    check(bool(input), "Truncated VICAR header.");
    const auto headerSize = number(vicar(std::string(lead, sizeof lead), "LBLSIZE"));
    check(headerSize >= sizeof lead && headerSize <= (firstRecord-1)*record && headerSize%record == 0,
          "Invalid VICAR header size.");
    input.seekg(0);
    std::string header(headerSize, '\0');
    input.read(header.data(), static_cast<std::streamsize>(header.size()));
    check(bool(input), "Truncated VICAR header.");
    check(vicar(header, "FORMAT") == "HALF" && vicar(header, "ORG") == "BSQ" &&
          vicar(header, "INTFMT") == "HIGH" && number(vicar(header, "NB")) == 1 &&
          number(vicar(header, "EOL")) == 0, "Unsupported VICAR pixel organization.");
    check(number(vicar(header, "NL")) == lines && number(vicar(header, "NS")) == samples &&
          number(vicar(header, "NBB")) == prefix && number(vicar(header, "RECSIZE")) == record &&
          number(vicar(header, "NLB")) == firstRecord-1-headerSize/record,
          "PDS and VICAR image layouts disagree.");

    CassiniImage result;
    result.imageMidTimeUtc = field(label, "IMAGE_MID_TIME");
    result.targetDescription = field(label, "TARGET_DESC");
    std::smatch gain;
    const auto gainText = field(label, "GAIN_MODE_ID");
    check(std::regex_match(gainText, gain, std::regex("([0-9.]+) ELECTRONS PER DN")),
          "Unsupported detector gain label.");
    result.electronsPerDn = std::stod(gain[1]);
    result.exposureMilliseconds = std::stod(field(label, "EXPOSURE_DURATION"));
    check(std::isfinite(result.electronsPerDn) && result.electronsPerDn > 0 &&
          std::isfinite(result.exposureMilliseconds) && result.exposureMilliseconds > 0,
          "Invalid exposure or gain.");
    result.dn.resize(lines, samples);
    input.seekg(static_cast<std::streamoff>((firstRecord-1)*record));
    std::vector<unsigned char> row(record);
    for (std::size_t y = 0; y < lines; ++y) {
        input.read(reinterpret_cast<char*>(row.data()), static_cast<std::streamsize>(record));
        check(bool(input), "Truncated IMG pixel record.");
        for (std::size_t x = 0; x < samples; ++x) {
            const unsigned value = (unsigned(row[prefix+2*x]) << 8) | row[prefix+2*x+1];
            const int signedValue = value < 32768 ? int(value) : int(value)-65536;
            result.dn(y, x) = signedValue;
        }
    }
    return result;
}
} // namespace fd::opnav::image
