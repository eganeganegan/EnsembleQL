#pragma once

#include "ensembleql/event.hpp"
#include "ensembleql/parser.hpp"
#include "ensembleql/planner.hpp"
#include "ensembleql/trajectory.hpp"

#include <string>
#include <vector>

namespace ensembleql {

class Engine {
public:
    std::vector<Event> query(Trajectory& trajectory, const std::string& query_text) const;
    std::vector<Event> execute(Trajectory& trajectory, const ExecutionPlan& plan) const;
};

} // namespace ensembleql
