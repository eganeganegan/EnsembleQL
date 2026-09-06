#pragma once

#include "ensembleql/ast.hpp"
#include "ensembleql/topology.hpp"

#include <memory>
#include <string>
#include <vector>

namespace ensembleql {

struct PlanNode {
    ast::Kind kind;
    std::shared_ptr<ast::Expr> expression;
    std::vector<std::shared_ptr<PlanNode>> inputs;
};

struct ExecutionPlan {
    std::shared_ptr<PlanNode> root;
    std::vector<std::string> required_selections;
    std::vector<std::string> required_observables;
    std::size_t unique_observables{};
};

struct PlanExplanation {
    std::vector<std::string> selections;
    std::vector<std::string> observables;
    std::vector<std::string> frame_predicates;
    std::vector<std::string> temporal_operations;
    std::string tree;
    bool streaming{true};
};

class Planner {
public:
    ExecutionPlan plan(const ast::Query& query, const Topology& topology) const;
    PlanExplanation explain(const ExecutionPlan& plan) const;
};

} // namespace ensembleql
