#pragma once

#include <string>

namespace ensembleql {

struct Distance {
    double nm{};
};

struct Duration {
    double ps{};
};

struct Angle {
    double degrees{};
};

struct Area {
    double nm2{};
};

Distance parse_distance(const std::string& text);
Duration parse_duration(const std::string& text);
Angle parse_angle(const std::string& text);
Area parse_area(const std::string& text);

} // namespace ensembleql
