#pragma once

#include "ensembleql/event.hpp"

#include <optional>
#include <vector>

namespace ensembleql::temporal {

bool before(const Event& a, const Event& b);
bool after(const Event& a, const Event& b);
bool overlaps(const Event& a, const Event& b);
std::vector<Event> join_before(const std::vector<Event>& left, const std::vector<Event>& right);
std::vector<Event> join_after(const std::vector<Event>& left, const std::vector<Event>& right);
std::vector<Event> join_overlaps(const std::vector<Event>& left, const std::vector<Event>& right);
std::vector<Event> followed_by(const std::vector<Event>& left,
                               const std::vector<Event>& right,
                               std::optional<double> within_ps = std::nullopt);

} // namespace ensembleql::temporal
