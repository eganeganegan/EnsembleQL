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
    std::size_t unique_observables{};
};

class Planner {
public:
    ExecutionPlan plan(const ast::Query& query, const Topology& topology) const;
};

} // namespace ensembleql
