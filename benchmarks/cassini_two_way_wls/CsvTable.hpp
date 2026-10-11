#pragma once
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <stdexcept>
#include <vector>

namespace cassini_benchmark {
// Reader for DeepNAV numeric benchmark exports, not a general CSV parser.
class CsvTable {
public:
  std::vector<std::string> header;
  std::vector<std::vector<std::string>> rows;
  explicit CsvTable(const std::filesystem::path &path) {
    std::ifstream input(path);
    if (!input)
      throw std::runtime_error("Cannot read " + path.string());
    std::string line;
    const auto split = [](const std::string &text) {
      if (text.find('"') != std::string::npos)
        throw std::runtime_error(
            "Quoted fields unsupported in numeric benchmark CSV");
      std::vector<std::string> out;
      std::size_t start = 0;
      for (;;) {
        auto end = text.find(',', start);
        out.push_back(text.substr(start, end - start));
        if (end == std::string::npos)
          break;
        start = end + 1;
      }
      return out;
    };
    if (!std::getline(input, line))
      throw std::runtime_error("Empty CSV");
    header = split(line);
    while (std::getline(input, line)) {
      auto row = split(line);
      if (row.size() != header.size())
        throw std::runtime_error("CSV row width mismatch");
      rows.push_back(std::move(row));
    }
  }
  std::size_t column(const std::string &name) const {
    for (std::size_t i = 0; i < header.size(); ++i)
      if (header[i] == name)
        return i;
    throw std::runtime_error("Missing CSV column " + name);
  }
  double number(const std::vector<std::string> &row,
                const std::string &name) const {
    return std::stod(row.at(column(name)));
  }
  void set(std::vector<std::string> &row, const std::string &name,
           double value) const {
    std::ostringstream out;
    out << std::setprecision(17) << value;
    row.at(column(name)) = out.str();
  }
  void write(const std::filesystem::path &path) const {
    std::ofstream out(path);
    out.exceptions(std::ios::badbit | std::ios::failbit);
    const auto record = [&](const auto &fields) {
      for (std::size_t i = 0; i < fields.size(); ++i)
        out << (i ? "," : "") << fields[i];
      out << '\n';
    };
    record(header);
    for (const auto &row : rows)
      record(row);
  }
};
} // namespace cassini_benchmark
