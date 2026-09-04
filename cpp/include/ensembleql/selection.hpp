#pragma once

#include "ensembleql/topology.hpp"

#include <cstddef>
#include <string>
#include <vector>

namespace ensembleql {

class Selection {
public:
    explicit Selection(std::string expression);
    const std::string& expression() const noexcept { return expression_; }
    std::vector<std::size_t> resolve(const Topology& topology) const;
private:
    std::string expression_;
};

} // namespace ensembleql
