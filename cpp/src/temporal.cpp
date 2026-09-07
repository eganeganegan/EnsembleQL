#include "ensembleql/temporal.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>
#include <utility>

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

bool touching(double first, double second) {
    const double tolerance = 1e-12 * std::max({1.0, std::abs(first), std::abs(second)});
    return std::abs(first - second) <= tolerance;
}

template <typename Predicate>
std::vector<Event> join(const std::vector<Event>& left, const std::vector<Event>& right,
                        const std::string& relation, Predicate predicate) {
    std::vector<Event> result;
    for (const auto& a : left) for (const auto& b : right) if (predicate(a, b)) result.push_back(combined(a, b, relation));
    return result;
}

std::vector<Event> earliest_right(const std::vector<Event>& left,
                                  const std::vector<Event>& right,
                                  const std::string& relation,
                                  std::optional<double> within_ps = std::nullopt) {
    std::vector<Event> result;
    for (const auto& a : left) {
        const Event* best = nullptr;
        for (const auto& b : right) {
            if (!before(a, b)) continue;
            const double gap = b.start_time - a.end_time;
            if (within_ps && gap > *within_ps) continue;
            if (!best || b.start_time < best->start_time) best = &b;
        }
        if (best) result.push_back(combined(a, *best, relation));
    }
    return result;
}
} // namespace

bool before(const Event& a, const Event& b) { return a.end_time <= b.start_time; }
bool after(const Event& a, const Event& b) { return a.start_time >= b.end_time; }
bool overlaps(const Event& a, const Event& b) { return a.start_time <= b.end_time && b.start_time <= a.end_time; }
bool during(const Event& a, const Event& b) {
    return b.start_time <= a.start_time && a.end_time <= b.end_time;
}
bool precedes(const Event& a, const Event& b) { return a.end_time < b.start_time; }
bool immediately_followed_by(const Event& a, const Event& b) {
    return before(a, b) && touching(a.end_time, b.start_time);
}
std::vector<Event> join_before(const std::vector<Event>& a, const std::vector<Event>& b) { return join(a, b, "BEFORE", before); }
std::vector<Event> join_after(const std::vector<Event>& a, const std::vector<Event>& b) { return join(a, b, "AFTER", after); }
std::vector<Event> join_overlaps(const std::vector<Event>& a, const std::vector<Event>& b) { return join(a, b, "OVERLAPS", overlaps); }
std::vector<Event> join_during(const std::vector<Event>& a, const std::vector<Event>& b) { return join(a, b, "DURING", during); }
std::vector<Event> join_precedes(const std::vector<Event>& a, const std::vector<Event>& b) { return join(a, b, "PRECEDES", precedes); }

std::vector<Event> followed_by(const std::vector<Event>& left, const std::vector<Event>& right,
                               std::optional<double> within_ps) {
    return earliest_right(left, right, "FOLLOWED_BY", within_ps);
}

std::vector<Event> join_immediately_followed_by(const std::vector<Event>& left,
                                                const std::vector<Event>& right) {
    return join(left, right, "IMMEDIATELY_FOLLOWED_BY", immediately_followed_by);
}

std::vector<Event> join_until(const std::vector<Event>& left,
                              const std::vector<Event>& right) {
    return earliest_right(left, right, "UNTIL");
}

std::vector<Event> repeats(const std::vector<Event>& events, std::size_t count,
                           std::optional<double> within_ps) {
    if (count == 0) throw std::invalid_argument("REPEATS count must be positive");
    std::vector<Event> result;
    if (events.size() < count) return result;
    for (std::size_t first = 0; first + count <= events.size(); ++first) {
        const Event& start = events[first];
        const Event& end = events[first + count - 1];
        const double span = end.end_time - start.start_time;
        if (within_ps && span > *within_ps) continue;
        Event repeated{start.start_time, end.end_time, "REPEATS", start.selections, {}};
        repeated.metadata["source_type"] = start.type;
        repeated.metadata["repeat_count"] = std::to_string(count);
        repeated.metadata["repeat_span_ps"] = std::to_string(span);
        result.push_back(std::move(repeated));
    }
    return result;
}

} // namespace ensembleql::temporal
