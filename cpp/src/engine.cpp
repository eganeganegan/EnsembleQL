#include "ensembleql/engine.hpp"

#include "ensembleql/observable.hpp"
#include "ensembleql/selection.hpp"
#include "ensembleql/temporal.hpp"

#include <algorithm>
#include <cmath>
#include <functional>
#include <map>
#include <memory>
#include <optional>
#include <sstream>
#include <stdexcept>

namespace ensembleql {
namespace {
using ValueCache = std::map<std::string, double>;
using ObservableRegistry = std::map<std::string, std::shared_ptr<Observable>>;
using Predicate = std::function<bool(const Frame&, ValueCache&)>;

std::string number(double value) {
    std::ostringstream output;
    output << value;
    return output.str();
}

std::string pair_signature(std::string name, std::string a, std::string b, const std::string& option = {}) {
    if (b < a) std::swap(a, b);
    return name + "(" + a + "," + b + option + ")";
}

std::string observable_signature(const ast::ExprPtr& expression) {
    switch (expression->kind) {
        case ast::Kind::Contact: {
            const auto& x = static_cast<const ast::ContactExpr&>(*expression);
            return pair_signature("CONTACT", x.a, x.b, ",cutoff=" + number(x.cutoff.nm) + "nm");
        }
        case ast::Kind::Distance: {
            const auto& x = static_cast<const ast::DistanceExpr&>(*expression);
            return pair_signature("DISTANCE", x.a, x.b);
        }
        case ast::Kind::ContactCount: {
            const auto& x = static_cast<const ast::ContactCountExpr&>(*expression);
            const std::string mode = x.mode == ContactMode::Residue ? ",mode=residue" : ",mode=atom";
            return pair_signature("CONTACT_COUNT", x.a, x.b,
                                  ",cutoff=" + number(x.cutoff.nm) + "nm" + mode);
        }
        case ast::Kind::Rg: {
            const auto& x = static_cast<const ast::RgExpr&>(*expression);
            return "RG(" + x.selection + (x.mass_weighted ? ",mass_weighted=true" : "") + ")";
        }
        case ast::Kind::HydrogenBond: {
            const auto& x = static_cast<const ast::HydrogenBondExpr&>(*expression);
            return "HBOND(" + x.donors + "," + x.acceptors +
                ",distance=" + number(x.distance_cutoff.nm) + "nm,angle=" +
                number(x.minimum_angle.degrees) + "deg)";
        }
        case ast::Kind::Dihedral: {
            const auto& x = static_cast<const ast::DihedralExpr&>(*expression);
            return "DIHEDRAL(" + x.selections[0] + "," + x.selections[1] + "," +
                   x.selections[2] + "," + x.selections[3] + ")";
        }
        case ast::Kind::Helix: {
            const auto& x = static_cast<const ast::HelixExpr&>(*expression);
            return "HELIX(" + x.selection + ",minimum_fraction=" + number(x.minimum_fraction) + ")";
        }
        case ast::Kind::Sasa: {
            const auto& x = static_cast<const ast::SasaExpr&>(*expression);
            return "SASA(" + x.selection + ",probe=" + number(x.probe.nm) +
                   "nm,points=" + std::to_string(x.points) + ")";
        }
        case ast::Kind::Rmsd: {
            const auto& x = static_cast<const ast::RmsdExpr&>(*expression);
            return "RMSD(" + x.selection + ",align=" + (x.align ? "true" : "false") + ")";
        }
        case ast::Kind::CoordinationNumber: {
            const auto& x = static_cast<const ast::CoordinationNumberExpr&>(*expression);
            return "COORDINATION_NUMBER(" + x.centers + "," + x.neighbors +
                   ",cutoff=" + number(x.cutoff.nm) + "nm)";
        }
        case ast::Kind::SaltBridge: {
            const auto& x = static_cast<const ast::SaltBridgeExpr&>(*expression);
            return pair_signature("SALT_BRIDGE", x.positive, x.negative,
                                  ",cutoff=" + number(x.cutoff.nm) + "nm");
        }
        case ast::Kind::AromaticStacking: {
            const auto& x = static_cast<const ast::AromaticStackingExpr&>(*expression);
            return pair_signature("AROMATIC_STACKING", x.first_ring, x.second_ring,
                ",distance=" + number(x.centroid_cutoff.nm) + "nm,angle=" +
                number(x.maximum_angle.degrees) + "deg");
        }
        case ast::Kind::SurfaceDistance: {
            const auto& x = static_cast<const ast::SurfaceDistanceExpr&>(*expression);
            return pair_signature("SURFACE_DISTANCE", x.molecule, x.surface);
        }
        case ast::Kind::Orientation: {
            const auto& x = static_cast<const ast::OrientationExpr&>(*expression);
            return "ORIENTATION(" + x.origin + "," + x.target + ",axis=" +
                   std::string(1, x.axis) + ")";
        }
        default:
            throw std::logic_error("Expression is not an observable");
    }
}

std::string signature(const ast::ExprPtr& expression) {
    switch (expression->kind) {
        case ast::Kind::Contact: {
            const auto& x = static_cast<const ast::ContactExpr&>(*expression);
            const std::string option = x.cutoff.nm == 0.45 ? "" : ",cutoff=" + number(x.cutoff.nm) + "nm";
            return "CONTACT(" + x.a + "," + x.b + option + ")";
        }
        case ast::Kind::Distance: { const auto& x = static_cast<const ast::DistanceExpr&>(*expression); return "DISTANCE(" + x.a + "," + x.b + ")"; }
        case ast::Kind::ContactCount: {
            const auto& x = static_cast<const ast::ContactCountExpr&>(*expression);
            std::string option = x.cutoff.nm == 0.45 ? "" : ",cutoff=" + number(x.cutoff.nm) + "nm";
            if (x.mode == ContactMode::Residue) option += ",mode=residue";
            return "CONTACT_COUNT(" + x.a + "," + x.b + option + ")";
        }
        case ast::Kind::Rg: {
            const auto& x = static_cast<const ast::RgExpr&>(*expression);
            return "RG(" + x.selection + (x.mass_weighted ? ",mass_weighted=true" : "") + ")";
        }
        case ast::Kind::HydrogenBond: case ast::Kind::Helix:
        case ast::Kind::SaltBridge: case ast::Kind::AromaticStacking:
            return observable_signature(expression);
        case ast::Kind::Dihedral: case ast::Kind::Sasa: case ast::Kind::Rmsd:
        case ast::Kind::CoordinationNumber: case ast::Kind::SurfaceDistance:
        case ast::Kind::Orientation:
            return observable_signature(expression);
        case ast::Kind::Comparison: { const auto& x = static_cast<const ast::ComparisonExpr&>(*expression); return signature(x.operand) + "#CMP" + std::to_string(static_cast<int>(x.op)) + ":" + std::to_string(x.threshold); }
        case ast::Kind::And: case ast::Kind::Or: { const auto& x = static_cast<const ast::BinaryExpr&>(*expression); return "(" + signature(x.left) + (expression->kind == ast::Kind::And ? "&" : "|") + signature(x.right) + ")"; }
        case ast::Kind::For: return signature(static_cast<const ast::ForExpr&>(*expression).operand);
        case ast::Kind::Repeats: return signature(static_cast<const ast::RepeatsExpr&>(*expression).operand);
        default: return "TEMPORAL";
    }
}

std::vector<std::string> expression_selections(const ast::ExprPtr& expression) {
    switch (expression->kind) {
        case ast::Kind::Contact: { const auto& x = static_cast<const ast::ContactExpr&>(*expression); return {x.a, x.b}; }
        case ast::Kind::Distance: { const auto& x = static_cast<const ast::DistanceExpr&>(*expression); return {x.a, x.b}; }
        case ast::Kind::ContactCount: { const auto& x = static_cast<const ast::ContactCountExpr&>(*expression); return {x.a, x.b}; }
        case ast::Kind::Rg: return {static_cast<const ast::RgExpr&>(*expression).selection};
        case ast::Kind::HydrogenBond: { const auto& x = static_cast<const ast::HydrogenBondExpr&>(*expression); return {x.donors, x.acceptors}; }
        case ast::Kind::Dihedral: { const auto& x = static_cast<const ast::DihedralExpr&>(*expression); return std::vector<std::string>(x.selections.begin(), x.selections.end()); }
        case ast::Kind::Helix: return {static_cast<const ast::HelixExpr&>(*expression).selection};
        case ast::Kind::Sasa: return {static_cast<const ast::SasaExpr&>(*expression).selection};
        case ast::Kind::Rmsd: return {static_cast<const ast::RmsdExpr&>(*expression).selection};
        case ast::Kind::CoordinationNumber: { const auto& x = static_cast<const ast::CoordinationNumberExpr&>(*expression); return {x.centers, x.neighbors}; }
        case ast::Kind::SaltBridge: { const auto& x = static_cast<const ast::SaltBridgeExpr&>(*expression); return {x.positive, x.negative}; }
        case ast::Kind::AromaticStacking: { const auto& x = static_cast<const ast::AromaticStackingExpr&>(*expression); return {x.first_ring, x.second_ring}; }
        case ast::Kind::SurfaceDistance: { const auto& x = static_cast<const ast::SurfaceDistanceExpr&>(*expression); return {x.molecule, x.surface}; }
        case ast::Kind::Orientation: { const auto& x = static_cast<const ast::OrientationExpr&>(*expression); return {x.origin, x.target}; }
        case ast::Kind::Comparison: return expression_selections(static_cast<const ast::ComparisonExpr&>(*expression).operand);
        case ast::Kind::And: case ast::Kind::Or: {
            const auto& x = static_cast<const ast::BinaryExpr&>(*expression);
            auto values = expression_selections(x.left); auto right = expression_selections(x.right);
            values.insert(values.end(), right.begin(), right.end()); return values;
        }
        default: return {};
    }
}

std::shared_ptr<Observable> make_observable(const ast::ExprPtr& expression, const Topology& topology) {
    auto resolve = [&](const std::string& value) { return Selection(value).resolve(topology); };
    switch (expression->kind) {
        case ast::Kind::Contact: { const auto& x = static_cast<const ast::ContactExpr&>(*expression); return std::make_shared<ContactObservable>(resolve(x.a), resolve(x.b), x.cutoff.nm); }
        case ast::Kind::Distance: { const auto& x = static_cast<const ast::DistanceExpr&>(*expression); return std::make_shared<DistanceObservable>(resolve(x.a), resolve(x.b)); }
        case ast::Kind::ContactCount: {
            const auto& x = static_cast<const ast::ContactCountExpr&>(*expression);
            return std::make_shared<ContactCountObservable>(
                resolve(x.a), resolve(x.b), x.cutoff.nm, x.mode, topology);
        }
        case ast::Kind::Rg: {
            const auto& x = static_cast<const ast::RgExpr&>(*expression);
            return std::make_shared<RgObservable>(resolve(x.selection), topology.bonds(),
                                                  x.mass_weighted ? topology.masses_da() : std::vector<double>{});
        }
        case ast::Kind::HydrogenBond: {
            const auto& x = static_cast<const ast::HydrogenBondExpr&>(*expression);
            return std::make_shared<HydrogenBondObservable>(
                resolve(x.donors), resolve(x.acceptors), topology,
                x.distance_cutoff.nm, x.minimum_angle.degrees);
        }
        case ast::Kind::Dihedral: {
            const auto& x = static_cast<const ast::DihedralExpr&>(*expression);
            return std::make_shared<DihedralObservable>(
                resolve(x.selections[0]).front(), resolve(x.selections[1]).front(),
                resolve(x.selections[2]).front(), resolve(x.selections[3]).front());
        }
        case ast::Kind::Helix: {
            const auto& x = static_cast<const ast::HelixExpr&>(*expression);
            return std::make_shared<HelixObservable>(resolve(x.selection), topology,
                                                     x.minimum_fraction);
        }
        case ast::Kind::Sasa: {
            const auto& x = static_cast<const ast::SasaExpr&>(*expression);
            return std::make_shared<SasaObservable>(resolve(x.selection), topology,
                                                    x.probe.nm, x.points);
        }
        case ast::Kind::Rmsd: {
            const auto& x = static_cast<const ast::RmsdExpr&>(*expression);
            return std::make_shared<RmsdObservable>(resolve(x.selection), x.align);
        }
        case ast::Kind::CoordinationNumber: {
            const auto& x = static_cast<const ast::CoordinationNumberExpr&>(*expression);
            return std::make_shared<CoordinationNumberObservable>(
                resolve(x.centers), resolve(x.neighbors), x.cutoff.nm);
        }
        case ast::Kind::SaltBridge: {
            const auto& x = static_cast<const ast::SaltBridgeExpr&>(*expression);
            return std::make_shared<SaltBridgeObservable>(
                resolve(x.positive), resolve(x.negative), x.cutoff.nm);
        }
        case ast::Kind::AromaticStacking: {
            const auto& x = static_cast<const ast::AromaticStackingExpr&>(*expression);
            return std::make_shared<AromaticStackingObservable>(
                resolve(x.first_ring), resolve(x.second_ring), x.centroid_cutoff.nm,
                x.maximum_angle.degrees);
        }
        case ast::Kind::SurfaceDistance: {
            const auto& x = static_cast<const ast::SurfaceDistanceExpr&>(*expression);
            return std::make_shared<SurfaceDistanceObservable>(resolve(x.molecule),
                                                               resolve(x.surface));
        }
        case ast::Kind::Orientation: {
            const auto& x = static_cast<const ast::OrientationExpr&>(*expression);
            return std::make_shared<OrientationObservable>(resolve(x.origin), resolve(x.target),
                                                           x.axis);
        }
        default: throw std::logic_error("Expression is not an observable");
    }
}

Predicate make_predicate(const ast::ExprPtr& expression, const Topology& topology,
                         ObservableRegistry& registry) {
    const auto bind_observable = [&](const ast::ExprPtr& observable_expression) {
        const std::string key = observable_signature(observable_expression);
        auto [iterator, inserted] = registry.emplace(key, nullptr);
        if (inserted) iterator->second = make_observable(observable_expression, topology);
        return std::pair{key, iterator->second};
    };
    if (expression->kind == ast::Kind::Contact ||
        expression->kind == ast::Kind::HydrogenBond || expression->kind == ast::Kind::Helix ||
        expression->kind == ast::Kind::SaltBridge ||
        expression->kind == ast::Kind::AromaticStacking) {
        auto [key, observable] = bind_observable(expression);
        return [key = std::move(key), observable = std::move(observable)](const Frame& frame, ValueCache& cache) {
            const auto found = cache.find(key);
            if (found != cache.end()) return found->second != 0.0;
            const double value = observable->evaluate(frame);
            cache.emplace(key, value);
            return value != 0.0;
        };
    }
    if (expression->kind == ast::Kind::Comparison) {
        const auto& comparison = static_cast<const ast::ComparisonExpr&>(*expression);
        auto [key, observable] = bind_observable(comparison.operand);
        const double threshold = comparison.threshold; const auto op = comparison.op;
        return [key = std::move(key), observable = std::move(observable), threshold, op](const Frame& frame, ValueCache& cache) {
            const auto found = cache.find(key);
            const double value = found != cache.end() ? found->second : observable->evaluate(frame);
            if (found == cache.end()) cache.emplace(key, value);
            switch (op) {
                case ast::ComparisonOp::Less: return value < threshold;
                case ast::ComparisonOp::LessEqual: return value <= threshold;
                case ast::ComparisonOp::Greater: return value > threshold;
                case ast::ComparisonOp::GreaterEqual: return value >= threshold;
                case ast::ComparisonOp::Equal: return value == threshold;
                case ast::ComparisonOp::NotEqual: return value != threshold;
            }
            return false;
        };
    }
    if (expression->kind == ast::Kind::And || expression->kind == ast::Kind::Or) {
        const auto& binary = static_cast<const ast::BinaryExpr&>(*expression);
        auto left = make_predicate(binary.left, topology, registry);
        auto right = make_predicate(binary.right, topology, registry);
        const bool is_and = expression->kind == ast::Kind::And;
        return [left = std::move(left), right = std::move(right), is_and](const Frame& frame, ValueCache& cache) {
            const bool left_value = left(frame, cache);
            const bool right_value = right(frame, cache);
            return is_and ? left_value && right_value : left_value || right_value;
        };
    }
    throw QueryError("Expression does not produce frame predicates");
}

void collect_sources(const ast::ExprPtr& expression, std::map<std::string, ast::ExprPtr>& sources) {
    if (expression->kind == ast::Kind::For) { collect_sources(static_cast<const ast::ForExpr&>(*expression).operand, sources); return; }
    if (expression->kind == ast::Kind::Repeats) { collect_sources(static_cast<const ast::RepeatsExpr&>(*expression).operand, sources); return; }
    if (expression->kind == ast::Kind::FollowedBy ||
        expression->kind == ast::Kind::ImmediatelyFollowedBy ||
        expression->kind == ast::Kind::Overlaps || expression->kind == ast::Kind::During ||
        expression->kind == ast::Kind::Until || expression->kind == ast::Kind::Before ||
        expression->kind == ast::Kind::Precedes || expression->kind == ast::Kind::After) {
        const auto& binary = static_cast<const ast::BinaryExpr&>(*expression);
        collect_sources(binary.left, sources); collect_sources(binary.right, sources); return;
    }
    sources.emplace(signature(expression), expression);
}

std::vector<Event> evaluate_events(const ast::ExprPtr& expression,
                                   const std::map<std::string, std::vector<Event>>& streams) {
    if (expression->kind == ast::Kind::For) {
        const auto& value = static_cast<const ast::ForExpr&>(*expression);
        return filter_for(evaluate_events(value.operand, streams), value.duration.ps);
    }
    if (expression->kind == ast::Kind::Repeats) {
        const auto& value = static_cast<const ast::RepeatsExpr&>(*expression);
        return temporal::repeats(evaluate_events(value.operand, streams), value.count,
                                 value.within ? std::optional<double>(value.within->ps) : std::nullopt);
    }
    if (expression->kind == ast::Kind::FollowedBy) {
        const auto& x = static_cast<const ast::FollowedByExpr&>(*expression);
        return temporal::followed_by(evaluate_events(x.left, streams), evaluate_events(x.right, streams), x.within ? std::optional<double>(x.within->ps) : std::nullopt);
    }
    if (expression->kind == ast::Kind::ImmediatelyFollowedBy ||
        expression->kind == ast::Kind::Overlaps || expression->kind == ast::Kind::During ||
        expression->kind == ast::Kind::Until || expression->kind == ast::Kind::Before ||
        expression->kind == ast::Kind::Precedes || expression->kind == ast::Kind::After) {
        const auto& x = static_cast<const ast::BinaryExpr&>(*expression);
        const auto left = evaluate_events(x.left, streams), right = evaluate_events(x.right, streams);
        if (expression->kind == ast::Kind::ImmediatelyFollowedBy) {
            return temporal::join_immediately_followed_by(left, right);
        }
        if (expression->kind == ast::Kind::Overlaps) return temporal::join_overlaps(left, right);
        if (expression->kind == ast::Kind::During) return temporal::join_during(left, right);
        if (expression->kind == ast::Kind::Until) return temporal::join_until(left, right);
        if (expression->kind == ast::Kind::Before) return temporal::join_before(left, right);
        if (expression->kind == ast::Kind::Precedes) return temporal::join_precedes(left, right);
        return temporal::join_after(left, right);
    }
    return streams.at(signature(expression));
}
} // namespace

std::vector<Event> Engine::query(Trajectory& trajectory, const std::string& query_text) const {
    return execute(trajectory, compile(trajectory.topology(), query_text));
}

PlanExplanation Engine::explain(const Topology& topology, const std::string& query_text) const {
    return Planner().explain(compile(topology, query_text));
}

const ExecutionPlan& Engine::compile(const Topology& topology, const std::string& query_text) const {
    if (cached_topology_identity_ != topology.cache_identity()) {
        plan_cache_.clear();
        cached_topology_identity_ = topology.cache_identity();
    }
    const auto found = plan_cache_.find(query_text);
    if (found != plan_cache_.end()) return found->second;
    auto [iterator, inserted] = plan_cache_.emplace(
        query_text, Planner().plan(Parser().parse(query_text), topology));
    (void)inserted;
    return iterator->second;
}

void Engine::clear_plan_cache() const noexcept {
    plan_cache_.clear();
    cached_topology_identity_ = 0;
}

std::vector<Event> Engine::execute(Trajectory& trajectory, const ExecutionPlan& plan) const {
    if (!plan.root || !plan.root->expression) throw std::invalid_argument("Execution plan has no root");
    std::map<std::string, ast::ExprPtr> sources;
    collect_sources(plan.root->expression, sources);
    struct State { Predicate predicate; EventExtractor extractor; std::vector<Event> events; };
    std::map<std::string, State> states;
    ObservableRegistry observables;
    for (const auto& [key, expression] : sources) {
        states.emplace(key, State{make_predicate(expression, trajectory.topology(), observables), EventExtractor(key, expression_selections(expression)), {}});
    }
    trajectory.reader().reset();
    Frame frame;
    std::optional<double> previous_time;
    while (trajectory.reader().next(frame)) {
        if (!std::isfinite(frame.time_ps)) throw std::runtime_error("Trajectory frame time must be finite");
        if (previous_time && frame.time_ps <= *previous_time) {
            throw std::runtime_error("Trajectory frame times must be strictly increasing (received " +
                                     number(frame.time_ps) + " ps after " + number(*previous_time) + " ps)");
        }
        previous_time = frame.time_ps;
        ValueCache values;
        for (auto& [key, state] : states) {
            (void)key;
            if (auto event = state.extractor.push(frame.time_ps, state.predicate(frame, values))) state.events.push_back(std::move(*event));
        }
    }
    std::map<std::string, std::vector<Event>> streams;
    for (auto& [key, state] : states) {
        if (auto event = state.extractor.finish()) state.events.push_back(std::move(*event));
        streams.emplace(key, std::move(state.events));
    }
    return evaluate_events(plan.root->expression, streams);
}

} // namespace ensembleql
