#include "ensembleql/event.hpp"
#include "ensembleql/geometry.hpp"
#include "ensembleql/temporal.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <iostream>
#include <vector>

using namespace ensembleql;

namespace {
std::size_t brute_force_contact_count(const Frame& frame,
                                      const std::vector<std::size_t>& left,
                                      const std::vector<std::size_t>& right,
                                      double cutoff_nm) {
    std::size_t count = 0;
    const double tolerance = 1e-12 * std::max(1.0, std::abs(cutoff_nm));
    for (const std::size_t i : left) {
        for (const std::size_t j : right) {
            if (distance(frame.coordinates[i], frame.coordinates[j]) <= cutoff_nm + tolerance) ++count;
        }
    }
    return count;
}
} // namespace

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
        const auto optimized_ms = std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - start).count();
        std::reverse(left.begin(), left.end());
        const auto boolean_start = std::chrono::steady_clock::now();
        const bool any = has_contact(frame, left, right, 0.45);
        const auto boolean_ms = std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - boolean_start).count();
        const auto brute_start = std::chrono::steady_clock::now();
        const std::size_t brute_count = brute_force_contact_count(frame, left, right, 0.45);
        const auto brute_ms = std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - brute_start).count();
        std::cout << "contacts atoms=" << atoms << " pairs=" << left.size() * right.size()
                  << " found=" << found.size() << " reference_found=" << brute_count
                  << " optimized_ms=" << optimized_ms << " brute_ms=" << brute_ms
                  << " speedup=" << brute_ms / optimized_ms
                  << " any=" << any << " boolean_ms=" << boolean_ms << '\n';
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
