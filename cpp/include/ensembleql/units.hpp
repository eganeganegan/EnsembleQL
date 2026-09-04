#pragma once

#include <string>

namespace ensembleql {

struct Distance {
    double nm{};
};

struct Duration {
    double ps{};
};

Distance parse_distance(const std::string& text);
Duration parse_duration(const std::string& text);

} // namespace ensembleql
