#include "ensembleql/event.hpp"

#include <utility>

namespace ensembleql {

EventExtractor::EventExtractor(std::string type, std::vector<std::string> selections)
    : type_(std::move(type)), selections_(std::move(selections)) {}

std::optional<Event> EventExtractor::push(double time_ps, bool predicate) {
    if (predicate) {
        if (!active_) { active_ = true; start_ = time_ps; }
        last_true_ = time_ps;
        return std::nullopt;
    }
    if (!active_) return std::nullopt;
    active_ = false;
    return Event{start_, last_true_, type_, selections_, {}};
}

std::optional<Event> EventExtractor::finish() {
    if (!active_) return std::nullopt;
    active_ = false;
    return Event{start_, last_true_, type_, selections_, {}};
}

std::vector<Event> filter_for(const std::vector<Event>& events, double minimum_duration_ps) {
    std::vector<Event> result;
    for (const auto& event : events) if (event.duration() >= minimum_duration_ps) result.push_back(event);
    return result;
}

} // namespace ensembleql
