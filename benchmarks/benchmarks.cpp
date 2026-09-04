#include "ensembleql/event.hpp"
#include "ensembleql/geometry.hpp"
#include "ensembleql/temporal.hpp"

#include <chrono>
#include <iostream>
#include <vector>

using namespace ensembleql;

int main() {
    for (const std::size_t atoms : {100U, 1000U, 5000U}) {
        Frame frame;
        frame.coordinates.resize(atoms);
        for (std::size_t i = 0; i < atoms; ++i) frame.coordinates[i] = {static_cast<double>(i) * 0.01, 0, 0};
        std::vector<std::size_t> left, right;
        for (std::size_t i = 0; i < atoms / 2; ++i) left.push_back(i);
        for (std::size_t i = atoms / 2; i < atoms; ++i) right.push_back(i);
        const auto start = std::chrono::steady_clock::now();
        const auto found = contacts(frame, left, right, 0.45);
        const auto elapsed = std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - start).count();
        std::cout << "contacts atoms=" << atoms << " pairs=" << left.size() * right.size()
                  << " found=" << found.size() << " ms=" << elapsed << '\n';
    }
    constexpr std::size_t count = 100000;
    EventExtractor extractor("synthetic");
    const auto extraction_start = std::chrono::steady_clock::now();
    std::vector<Event> events;
    for (std::size_t i = 0; i < count; ++i) if (auto event = extractor.push(static_cast<double>(i), i % 10 < 5)) events.push_back(*event);
    if (auto event = extractor.finish()) events.push_back(*event);
    const auto extraction_ms = std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - extraction_start).count();
    const auto join_start = std::chrono::steady_clock::now();
    const auto joined = temporal::followed_by(events, events, 10.0);
    const auto join_ms = std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - join_start).count();
    std::cout << "event_extraction frames=" << count << " events=" << events.size() << " ms=" << extraction_ms << '\n';
    std::cout << "temporal_join events=" << events.size() << " matches=" << joined.size() << " ms=" << join_ms << '\n';
}
