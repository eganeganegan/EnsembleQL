#include "ensembleql/temporal.hpp"

#include <algorithm>
#include <string>

namespace ensembleql::temporal {
namespace {
Event combined(const Event& a, const Event& b, const std::string& relation) {
    Event result{std::min(a.start_time, b.start_time), std::max(a.end_time, b.end_time), relation, {}, {}};
    result.metadata["left_type"] = a.type;
    result.metadata["right_type"] = b.type;
    result.metadata["transition_gap_ps"] = std::to_string(std::max(0.0, b.start_time - a.end_time));
    result.selections = a.selections;
    for (const auto& selection : b.selections) if (std::find(result.selections.begin(), result.selections.end(), selection) == result.selections.end()) result.selections.push_back(selection);
    return result;
}
template <typename Predicate>
std::vector<Event> join(const std::vector<Event>& left, const std::vector<Event>& right,
                        const std::string& relation, Predicate predicate) {
    std::vector<Event> result;
    for (const auto& a : left) for (const auto& b : right) if (predicate(a, b)) result.push_back(combined(a, b, relation));
    return result;
}
} // namespace

bool before(const Event& a, const Event& b) { return a.end_time <= b.start_time; }
bool after(const Event& a, const Event& b) { return a.start_time >= b.end_time; }
bool overlaps(const Event& a, const Event& b) { return a.start_time <= b.end_time && b.start_time <= a.end_time; }
std::vector<Event> join_before(const std::vector<Event>& a, const std::vector<Event>& b) { return join(a, b, "BEFORE", before); }
std::vector<Event> join_after(const std::vector<Event>& a, const std::vector<Event>& b) { return join(a, b, "AFTER", after); }
std::vector<Event> join_overlaps(const std::vector<Event>& a, const std::vector<Event>& b) { return join(a, b, "OVERLAPS", overlaps); }

std::vector<Event> followed_by(const std::vector<Event>& left, const std::vector<Event>& right,
                               std::optional<double> within_ps) {
    std::vector<Event> result;
    for (const auto& a : left) {
        const Event* best = nullptr;
        for (const auto& b : right) {
            if (!before(a, b)) continue;
            const double gap = b.start_time - a.end_time;
            if (within_ps && gap > *within_ps) continue; // WITHIN is inclusive
            if (!best || b.start_time < best->start_time) best = &b;
        }
        if (best) result.push_back(combined(a, *best, "FOLLOWED_BY"));
    }
    return result;
}

} // namespace ensembleql::temporal
