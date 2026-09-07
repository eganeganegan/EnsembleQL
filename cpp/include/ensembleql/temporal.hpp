#pragma once

#include "ensembleql/event.hpp"

#include <optional>
#include <vector>

namespace ensembleql::temporal {

bool before(const Event& a, const Event& b);
bool after(const Event& a, const Event& b);
bool overlaps(const Event& a, const Event& b);
bool during(const Event& a, const Event& b);
bool precedes(const Event& a, const Event& b);
bool immediately_followed_by(const Event& a, const Event& b);
std::vector<Event> join_before(const std::vector<Event>& left, const std::vector<Event>& right);
std::vector<Event> join_after(const std::vector<Event>& left, const std::vector<Event>& right);
std::vector<Event> join_overlaps(const std::vector<Event>& left, const std::vector<Event>& right);
std::vector<Event> join_during(const std::vector<Event>& left, const std::vector<Event>& right);
std::vector<Event> join_precedes(const std::vector<Event>& left, const std::vector<Event>& right);
std::vector<Event> followed_by(const std::vector<Event>& left,
                               const std::vector<Event>& right,
                               std::optional<double> within_ps = std::nullopt);
std::vector<Event> join_immediately_followed_by(const std::vector<Event>& left,
                                                const std::vector<Event>& right);
std::vector<Event> join_until(const std::vector<Event>& left,
                              const std::vector<Event>& right);
std::vector<Event> repeats(const std::vector<Event>& events,
                           std::size_t count,
                           std::optional<double> within_ps = std::nullopt);

} // namespace ensembleql::temporal
