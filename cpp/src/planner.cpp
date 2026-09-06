#include "ensembleql/planner.hpp"
#include "ensembleql/parser.hpp"
#include "ensembleql/selection.hpp"

#include <cmath>
#include <set>
#include <sstream>
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

std::string number(double value) {
    std::ostringstream output;
    output << value;
    return output.str();
}

std::string comparison_name(ast::ComparisonOp op) {
    switch (op) {
        case ast::ComparisonOp::Less: return "<";
        case ast::ComparisonOp::LessEqual: return "<=";
        case ast::ComparisonOp::Greater: return ">";
        case ast::ComparisonOp::GreaterEqual: return ">=";
        case ast::ComparisonOp::Equal: return "==";
        case ast::ComparisonOp::NotEqual: return "!=";
    }
    throw std::logic_error("Unknown comparison operator");
}

std::string render(const ast::ExprPtr& expression) {
    switch (expression->kind) {
        case ast::Kind::Contact: {
            const auto& value = static_cast<const ast::ContactExpr&>(*expression);
            return "CONTACT(" + value.a + ", " + value.b + ", cutoff=" + number(value.cutoff.nm) + "nm)";
        }
        case ast::Kind::Distance: {
            const auto& value = static_cast<const ast::DistanceExpr&>(*expression);
            return "DISTANCE(" + value.a + ", " + value.b + ")";
        }
        case ast::Kind::ContactCount: {
            const auto& value = static_cast<const ast::ContactCountExpr&>(*expression);
            return "CONTACT_COUNT(" + value.a + ", " + value.b + ", cutoff=" +
                   number(value.cutoff.nm) + "nm, mode=" +
                   (value.mode == ContactMode::Residue ? "residue" : "atom") + ")";
        }
        case ast::Kind::Rg: {
            const auto& value = static_cast<const ast::RgExpr&>(*expression);
            return "RG(" + value.selection +
                   (value.mass_weighted ? ", mass_weighted=true" : "") + ")";
        }
        case ast::Kind::Comparison: {
            const auto& value = static_cast<const ast::ComparisonExpr&>(*expression);
            return render(value.operand) + " " + comparison_name(value.op) + " " + number(value.threshold) +
                   (value.distance_threshold ? "nm" : "");
        }
        case ast::Kind::For: {
            const auto& value = static_cast<const ast::ForExpr&>(*expression);
            return render(value.operand) + " FOR >= " + number(value.duration.ps) + "ps";
        }
        case ast::Kind::FollowedBy: {
            const auto& value = static_cast<const ast::FollowedByExpr&>(*expression);
            return "(" + render(value.left) + " FOLLOWED_BY " + render(value.right) +
                   (value.within ? " WITHIN " + number(value.within->ps) + "ps" : "") + ")";
        }
        case ast::Kind::Overlaps: case ast::Kind::Before: case ast::Kind::After:
        case ast::Kind::And: case ast::Kind::Or: {
            const auto& value = static_cast<const ast::BinaryExpr&>(*expression);
            const std::string operation = expression->kind == ast::Kind::Overlaps ? "OVERLAPS" :
                expression->kind == ast::Kind::Before ? "BEFORE" :
                expression->kind == ast::Kind::After ? "AFTER" :
                expression->kind == ast::Kind::And ? "AND" : "OR";
            return "(" + render(value.left) + " " + operation + " " + render(value.right) + ")";
        }
    }
    throw std::logic_error("Unknown AST node");
}

std::string kind_name(ast::Kind kind) {
    switch (kind) {
        case ast::Kind::Contact: return "CONTACT";
        case ast::Kind::Distance: return "DISTANCE";
        case ast::Kind::ContactCount: return "CONTACT_COUNT";
        case ast::Kind::Rg: return "RG";
        case ast::Kind::Comparison: return "COMPARISON";
        case ast::Kind::For: return "FOR";
        case ast::Kind::FollowedBy: return "FOLLOWED_BY";
        case ast::Kind::Overlaps: return "OVERLAPS";
        case ast::Kind::Before: return "BEFORE";
        case ast::Kind::After: return "AFTER";
        case ast::Kind::And: return "AND";
        case ast::Kind::Or: return "OR";
    }
    throw std::logic_error("Unknown AST node");
}

void render_tree(const std::shared_ptr<PlanNode>& node, std::size_t depth, std::ostringstream& output) {
    output << std::string(depth * 2, ' ') << kind_name(node->kind) << ": " << render(node->expression) << '\n';
    for (const auto& input : node->inputs) render_tree(input, depth + 1, output);
}

void collect_explanation(const ast::ExprPtr& expression, std::vector<std::string>& predicates,
                         std::vector<std::string>& temporal) {
    if (expression->kind == ast::Kind::For) {
        const auto& value = static_cast<const ast::ForExpr&>(*expression);
        collect_explanation(value.operand, predicates, temporal);
        temporal.push_back("FOR >= " + number(value.duration.ps) + "ps");
        return;
    }
    if (expression->kind == ast::Kind::FollowedBy || expression->kind == ast::Kind::Overlaps ||
        expression->kind == ast::Kind::Before || expression->kind == ast::Kind::After) {
        const auto& value = static_cast<const ast::BinaryExpr&>(*expression);
        collect_explanation(value.left, predicates, temporal);
        collect_explanation(value.right, predicates, temporal);
        if (expression->kind == ast::Kind::FollowedBy) {
            const auto& followed = static_cast<const ast::FollowedByExpr&>(*expression);
            temporal.push_back("FOLLOWED_BY" + (followed.within ? " WITHIN " + number(followed.within->ps) + "ps" : ""));
        } else {
            temporal.push_back(kind_name(expression->kind));
        }
        return;
    }
    predicates.push_back(render(expression));
}

std::shared_ptr<PlanNode> build(const ast::ExprPtr& expression, const Topology& topology,
                                std::set<std::string>& selections, std::set<std::string>& observables) {
    if (!expression) throw std::invalid_argument("Cannot plan an empty query");
    auto node = std::make_shared<PlanNode>();
    node->kind = expression->kind;
    node->expression = expression;
    auto require = [&](const std::string& value) {
        auto resolved = Selection(value).resolve(topology);
        selections.insert(value);
        return resolved;
    };
    switch (expression->kind) {
        case ast::Kind::Contact: {
            const auto& value = static_cast<const ast::ContactExpr&>(*expression);
            require(value.a); require(value.b); observables.insert(pair_key("CONTACT", value.a, value.b, ",cutoff=" + number(value.cutoff.nm) + "nm")); break;
        }
        case ast::Kind::Distance: {
            const auto& value = static_cast<const ast::DistanceExpr&>(*expression);
            require(value.a); require(value.b); observables.insert(pair_key("DISTANCE", value.a, value.b)); break;
        }
        case ast::Kind::ContactCount: {
            const auto& value = static_cast<const ast::ContactCountExpr&>(*expression);
            const std::string mode = value.mode == ContactMode::Residue ? ",mode=residue" : ",mode=atom";
            require(value.a); require(value.b);
            observables.insert(pair_key("CONTACT_COUNT", value.a, value.b,
                                        ",cutoff=" + number(value.cutoff.nm) + "nm" + mode));
            break;
        }
        case ast::Kind::Rg: {
            const auto& value = static_cast<const ast::RgExpr&>(*expression);
            const auto resolved = require(value.selection);
            if (value.mass_weighted) {
                for (const std::size_t index : resolved) {
                    const double mass = topology.atoms()[index].mass_da;
                    if (!std::isfinite(mass) || mass <= 0.0) {
                        throw QueryError("Mass-weighted RG requires a known positive mass for atom " +
                                         std::to_string(index) + " (element '" +
                                         topology.atoms()[index].element + "')");
                    }
                }
            }
            observables.insert("RG(" + value.selection +
                               (value.mass_weighted ? ",mass_weighted=true" : "") + ")");
            break;
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
    return {std::move(root), {selections.begin(), selections.end()},
            {observables.begin(), observables.end()}, observables.size()};
}

PlanExplanation Planner::explain(const ExecutionPlan& plan) const {
    if (!plan.root || !plan.root->expression) throw std::invalid_argument("Execution plan has no root");
    PlanExplanation explanation;
    explanation.selections = plan.required_selections;
    explanation.observables = plan.required_observables;
    collect_explanation(plan.root->expression, explanation.frame_predicates, explanation.temporal_operations);
    std::ostringstream tree;
    render_tree(plan.root, 0, tree);
    explanation.tree = tree.str();
    return explanation;
}

} // namespace ensembleql
