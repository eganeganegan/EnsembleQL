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
std::string lowercase(std::string value) {
    std::transform(value.begin(), value.end(), value.begin(), [](unsigned char character) {
        return static_cast<char>(std::tolower(character));
    });
    return value;
}
} // namespace
Distance parse_distance(const std::string& text) {
    auto [value, unit] = split(text);
    const std::string normalized = lowercase(unit);
    if (normalized == "nm") return {value};
    if (normalized == "a" || normalized == "angstrom") return {value * 0.1};
    if (normalized == "fs" || normalized == "ps" || normalized == "ns" || normalized == "us") throw std::runtime_error("Distance expected, received time unit '" + unit + "'");
    throw std::runtime_error("Unknown distance unit '" + unit + "'");
}
Duration parse_duration(const std::string& text) {
    auto [value, unit] = split(text);
    const std::string normalized = lowercase(unit);
    if (normalized == "fs") return {value * 0.001};
    if (normalized == "ps") return {value};
    if (normalized == "ns") return {value * 1000.0};
    if (normalized == "us") return {value * 1000000.0};
    if (normalized == "nm" || normalized == "a" || normalized == "angstrom") throw std::runtime_error("Time expected, received distance unit '" + unit + "'");
    throw std::runtime_error("Unknown time unit '" + unit + "'");
}
} // namespace ensembleql
