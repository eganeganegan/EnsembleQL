#pragma once

#include <map>
#include <optional>
#include <string>
#include <vector>

namespace ensembleql {

struct Event {
    double start_time{};
    double end_time{};
    std::string type;
    std::vector<std::string> selections;
    std::map<std::string, std::string> metadata;
    double duration() const noexcept { return end_time - start_time; }
};

class EventExtractor {
public:
    explicit EventExtractor(std::string type, std::vector<std::string> selections = {});
    std::optional<Event> push(double time_ps, bool predicate);
    std::optional<Event> finish();
private:
    std::string type_;
    std::vector<std::string> selections_;
    bool active_{};
    double start_{};
    double last_true_{};
};

std::vector<Event> filter_for(const std::vector<Event>& events, double minimum_duration_ps);

} // namespace ensembleql
