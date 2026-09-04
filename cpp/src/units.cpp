#include "ensembleql/units.hpp"

#include <algorithm>
#include <cctype>
#include <regex>
#include <stdexcept>

namespace ensembleql {
namespace {
std::pair<double, std::string> split(const std::string& text) {
    const std::regex pattern(R"(^\s*([+-]?[0-9]+(?:\.[0-9]+)?)\s*([A-Za-z]+)\s*$)");
    std::smatch match;
    if (!std::regex_match(text, match, pattern)) throw std::runtime_error("Expected a number with an explicit unit, received '" + text + "'");
    return {std::stod(match[1].str()), match[2].str()};
}
} // namespace
Distance parse_distance(const std::string& text) {
    auto [value, unit] = split(text);
    if (unit == "nm") return {value};
    if (unit == "A" || unit == "angstrom" || unit == "Angstrom") return {value * 0.1};
    if (unit == "fs" || unit == "ps" || unit == "ns" || unit == "us") throw std::runtime_error("Distance expected, received time unit '" + unit + "'");
    throw std::runtime_error("Unknown distance unit '" + unit + "'");
}
Duration parse_duration(const std::string& text) {
    auto [value, unit] = split(text);
    if (unit == "fs") return {value * 0.001};
    if (unit == "ps") return {value};
    if (unit == "ns") return {value * 1000.0};
    if (unit == "us") return {value * 1000000.0};
    if (unit == "nm" || unit == "A" || unit == "angstrom") throw std::runtime_error("Time expected, received distance unit '" + unit + "'");
    throw std::runtime_error("Unknown time unit '" + unit + "'");
}
} // namespace ensembleql
