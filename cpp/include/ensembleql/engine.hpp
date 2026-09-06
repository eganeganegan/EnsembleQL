#pragma once

#include "ensembleql/event.hpp"
#include "ensembleql/parser.hpp"
#include "ensembleql/planner.hpp"
#include "ensembleql/trajectory.hpp"

#include <cstddef>
#include <cstdint>
#include <string>
#include <unordered_map>
#include <vector>

namespace ensembleql {

class Engine {
public:
    std::vector<Event> query(Trajectory& trajectory, const std::string& query_text) const;
    std::vector<Event> execute(Trajectory& trajectory, const ExecutionPlan& plan) const;
    PlanExplanation explain(const Topology& topology, const std::string& query_text) const;
    const ExecutionPlan& compile(const Topology& topology, const std::string& query_text) const;
    std::size_t cached_plan_count() const noexcept { return plan_cache_.size(); }
    void clear_plan_cache() const noexcept;
private:
    mutable std::uint64_t cached_topology_identity_{};
    mutable std::unordered_map<std::string, ExecutionPlan> plan_cache_;
};

} // namespace ensembleql
