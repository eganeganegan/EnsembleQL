#include "ensembleql/planner.hpp"
#include "ensembleql/parser.hpp"
#include "ensembleql/selection.hpp"

#include <set>
#include <stdexcept>
#include <utility>

namespace ensembleql {
namespace {
std::string pair_key(std::string name, std::string a, std::string b, const std::string& option = {}) {
    if (b < a) std::swap(a, b);
    return name + "(" + a + "," + b + option + ")";
}

bool frame_predicate(const ast::ExprPtr& expression) {
    if (expression->kind == ast::Kind::Contact || expression->kind == ast::Kind::Comparison) return true;
    if (expression->kind == ast::Kind::And || expression->kind == ast::Kind::Or) {
        const auto& value = static_cast<const ast::BinaryExpr&>(*expression);
        return frame_predicate(value.left) && frame_predicate(value.right);
    }
    return false;
}

bool event_expression(const ast::ExprPtr& expression) {
    if (frame_predicate(expression)) return true;
    if (expression->kind == ast::Kind::For) {
        return event_expression(static_cast<const ast::ForExpr&>(*expression).operand);
    }
    if (expression->kind == ast::Kind::FollowedBy || expression->kind == ast::Kind::Overlaps ||
        expression->kind == ast::Kind::Before || expression->kind == ast::Kind::After) {
        const auto& value = static_cast<const ast::BinaryExpr&>(*expression);
        return event_expression(value.left) && event_expression(value.right);
    }
    return false;
}

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
            require(value.a); require(value.b); observables.insert(pair_key("CONTACT", value.a, value.b, "," + std::to_string(value.cutoff.nm))); break;
        }
        case ast::Kind::Distance: {
            const auto& value = static_cast<const ast::DistanceExpr&>(*expression);
            require(value.a); require(value.b); observables.insert(pair_key("DISTANCE", value.a, value.b)); break;
        }
        case ast::Kind::ContactCount: {
            const auto& value = static_cast<const ast::ContactCountExpr&>(*expression);
            require(value.a); require(value.b); observables.insert(pair_key("CONTACT_COUNT", value.a, value.b, "," + std::to_string(value.cutoff.nm))); break;
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
            if (!event_expression(value.operand)) throw QueryError("FOR requires an event-producing expression");
            node->inputs.push_back(build(value.operand, topology, selections, observables)); break;
        }
        case ast::Kind::FollowedBy: case ast::Kind::Overlaps: case ast::Kind::Before:
        case ast::Kind::After: {
            const auto& value = static_cast<const ast::BinaryExpr&>(*expression);
            if (!event_expression(value.left) || !event_expression(value.right)) {
                const std::string name = expression->kind == ast::Kind::FollowedBy ? "FOLLOWED_BY" :
                    expression->kind == ast::Kind::Overlaps ? "OVERLAPS" :
                    expression->kind == ast::Kind::Before ? "BEFORE" : "AFTER";
                throw QueryError(name + " requires event-producing expressions");
            }
            node->inputs.push_back(build(value.left, topology, selections, observables));
            node->inputs.push_back(build(value.right, topology, selections, observables)); break;
        }
        case ast::Kind::And: case ast::Kind::Or: {
            const auto& value = static_cast<const ast::BinaryExpr&>(*expression);
            if (!frame_predicate(value.left) || !frame_predicate(value.right)) {
                throw QueryError(std::string(expression->kind == ast::Kind::And ? "AND" : "OR") +
                                 " requires frame-predicate expressions");
            }
            node->inputs.push_back(build(value.left, topology, selections, observables));
            node->inputs.push_back(build(value.right, topology, selections, observables)); break;
        }
    }
    return node;
}
} // namespace

ExecutionPlan Planner::plan(const ast::Query& query, const Topology& topology) const {
    if (!query.root || !event_expression(query.root)) {
        throw QueryError("Query root must produce events; numeric observables require a comparison");
    }
    std::set<std::string> selections, observables;
    auto root = build(query.root, topology, selections, observables);
    return {std::move(root), {selections.begin(), selections.end()}, observables.size()};
}

} // namespace ensembleql
