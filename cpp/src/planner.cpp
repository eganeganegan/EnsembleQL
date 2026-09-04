#include "ensembleql/planner.hpp"
#include "ensembleql/selection.hpp"

#include <set>
#include <stdexcept>

namespace ensembleql {
namespace {
std::shared_ptr<PlanNode> build(const ast::ExprPtr& expression, const Topology& topology,
                                std::set<std::string>& selections, std::set<std::string>& observables) {
    if (!expression) throw std::invalid_argument("Cannot plan an empty query");
    auto node = std::make_shared<PlanNode>();
    node->kind = expression->kind;
    node->expression = expression;
    auto require = [&](const std::string& value) {
        Selection(value).resolve(topology);
        selections.insert(value);
    };
    switch (expression->kind) {
        case ast::Kind::Contact: {
            const auto& value = static_cast<const ast::ContactExpr&>(*expression);
            require(value.a); require(value.b); observables.insert("CONTACT(" + value.a + "," + value.b + ")"); break;
        }
        case ast::Kind::Distance: {
            const auto& value = static_cast<const ast::DistanceExpr&>(*expression);
            require(value.a); require(value.b); observables.insert("DISTANCE(" + value.a + "," + value.b + ")"); break;
        }
        case ast::Kind::ContactCount: {
            const auto& value = static_cast<const ast::ContactCountExpr&>(*expression);
            require(value.a); require(value.b); observables.insert("CONTACT_COUNT(" + value.a + "," + value.b + ")"); break;
        }
        case ast::Kind::Rg: {
            const auto& value = static_cast<const ast::RgExpr&>(*expression);
            require(value.selection); observables.insert("RG(" + value.selection + ")"); break;
        }
        case ast::Kind::Comparison: {
            const auto& value = static_cast<const ast::ComparisonExpr&>(*expression);
            node->inputs.push_back(build(value.operand, topology, selections, observables)); break;
        }
        case ast::Kind::For: {
            const auto& value = static_cast<const ast::ForExpr&>(*expression);
            node->inputs.push_back(build(value.operand, topology, selections, observables)); break;
        }
        case ast::Kind::FollowedBy: case ast::Kind::Overlaps: case ast::Kind::Before:
        case ast::Kind::After: case ast::Kind::And: case ast::Kind::Or: {
            const auto& value = static_cast<const ast::BinaryExpr&>(*expression);
            node->inputs.push_back(build(value.left, topology, selections, observables));
            node->inputs.push_back(build(value.right, topology, selections, observables)); break;
        }
    }
    return node;
}
} // namespace

ExecutionPlan Planner::plan(const ast::Query& query, const Topology& topology) const {
    std::set<std::string> selections, observables;
    auto root = build(query.root, topology, selections, observables);
    return {std::move(root), {selections.begin(), selections.end()}, observables.size()};
}

} // namespace ensembleql
