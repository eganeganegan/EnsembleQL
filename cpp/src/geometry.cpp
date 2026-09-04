#include "ensembleql/geometry.hpp"

#include <cmath>
#include <limits>
#include <stdexcept>

namespace ensembleql {
namespace {
void validate(const Frame& frame, const std::vector<std::size_t>& indices) {
    if (frame.box_nm) throw std::runtime_error("Periodic boundary conditions are not yet supported");
    if (indices.empty()) throw std::invalid_argument("Geometry selection is empty");
    for (const auto index : indices) if (index >= frame.coordinates.size()) throw std::out_of_range("Atom index exceeds frame coordinates");
}
} // namespace

double distance(const Vec3& a, const Vec3& b) {
    const double dx = a[0] - b[0], dy = a[1] - b[1], dz = a[2] - b[2];
    return std::sqrt(dx * dx + dy * dy + dz * dz);
}

double minimum_distance(const Frame& frame, const std::vector<std::size_t>& a, const std::vector<std::size_t>& b) {
    validate(frame, a); validate(frame, b);
    double result = std::numeric_limits<double>::infinity();
    for (const auto i : a) for (const auto j : b) if (i != j) result = std::min(result, distance(frame.coordinates[i], frame.coordinates[j]));
    if (!std::isfinite(result)) throw std::invalid_argument("Selections contain no distinct atom pair");
    return result;
}

std::vector<Contact> contacts(const Frame& frame, const std::vector<std::size_t>& a,
                              const std::vector<std::size_t>& b, double cutoff_nm) {
    validate(frame, a); validate(frame, b);
    if (cutoff_nm < 0.0) throw std::invalid_argument("Contact cutoff must be non-negative");
    std::vector<Contact> result;
    for (const auto i : a) for (const auto j : b) {
        if (i == j) continue;
        const double d = distance(frame.coordinates[i], frame.coordinates[j]);
        if (d <= cutoff_nm) result.push_back({i, j, d}); // cutoff is inclusive
    }
    return result;
}

double radius_of_gyration(const Frame& frame, const std::vector<std::size_t>& selection) {
    validate(frame, selection);
    Vec3 center{};
    for (const auto i : selection) for (std::size_t d = 0; d < 3; ++d) center[d] += frame.coordinates[i][d];
    for (double& value : center) value /= static_cast<double>(selection.size());
    double squared = 0.0;
    for (const auto i : selection) {
        const double d = distance(frame.coordinates[i], center);
        squared += d * d;
    }
    return std::sqrt(squared / static_cast<double>(selection.size()));
}

} // namespace ensembleql
