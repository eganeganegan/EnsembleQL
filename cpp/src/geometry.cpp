#include "ensembleql/geometry.hpp"

#include <array>
#include <cmath>
#include <limits>
#include <optional>
#include <queue>
#include <stdexcept>

namespace ensembleql {
namespace {
void validate_indices(const Frame& frame, const std::vector<std::size_t>& indices) {
    if (indices.empty()) throw std::invalid_argument("Geometry selection is empty");
    for (const auto index : indices) if (index >= frame.coordinates.size()) throw std::out_of_range("Atom index exceeds frame coordinates");
}

void validate_box(const Vec3& box_nm) {
    for (const double length : box_nm) {
        if (!std::isfinite(length) || length <= 0.0) {
            throw std::invalid_argument("Orthorhombic box lengths must be finite and positive");
        }
    }
}

Vec3 cross(const Vec3& a, const Vec3& b) {
    return {a[1] * b[2] - a[2] * b[1],
            a[2] * b[0] - a[0] * b[2],
            a[0] * b[1] - a[1] * b[0]};
}

double dot(const Vec3& a, const Vec3& b) {
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

double squared_norm(const Vec3& value) {
    return dot(value, value);
}

void validate_cell(const PeriodicCell& cell) {
    for (const auto& vector : cell.vectors_nm) {
        for (const double value : vector) {
            if (!std::isfinite(value)) throw std::invalid_argument("Periodic cell vectors must be finite");
        }
    }
    const double determinant = dot(cell.vectors_nm[0], cross(cell.vectors_nm[1], cell.vectors_nm[2]));
    if (!std::isfinite(determinant) || std::abs(determinant) <= 1e-15) {
        throw std::invalid_argument("Periodic cell matrix must be finite and invertible");
    }
}

PeriodicCell orthorhombic_cell(const Vec3& box_nm) {
    validate_box(box_nm);
    PeriodicCell cell;
    cell.vectors_nm = {Vec3{box_nm[0], 0.0, 0.0},
                       Vec3{0.0, box_nm[1], 0.0},
                       Vec3{0.0, 0.0, box_nm[2]}};
    return cell;
}

std::optional<PeriodicCell> frame_cell(const Frame& frame) {
    if (frame.cell_nm) {
        validate_cell(*frame.cell_nm);
        return frame.cell_nm;
    }
    if (frame.box_nm) return orthorhombic_cell(*frame.box_nm);
    return std::nullopt;
}

double frame_distance(const std::optional<PeriodicCell>& cell, const Vec3& a, const Vec3& b) {
    return cell ? minimum_image_distance_cell(a, b, *cell) : distance(a, b);
}
} // namespace

double distance(const Vec3& a, const Vec3& b) {
    const double dx = a[0] - b[0], dy = a[1] - b[1], dz = a[2] - b[2];
    return std::sqrt(dx * dx + dy * dy + dz * dz);
}

double minimum_image_distance(const Vec3& a, const Vec3& b, const Vec3& box_nm) {
    return minimum_image_distance_cell(a, b, orthorhombic_cell(box_nm));
}

Vec3 minimum_image_displacement(const Vec3& displacement, const PeriodicCell& cell) {
    validate_cell(cell);
    const auto& a = cell.vectors_nm[0];
    const auto& b = cell.vectors_nm[1];
    const auto& c = cell.vectors_nm[2];
    const double determinant = dot(a, cross(b, c));
    const std::array<double, 3> fractional{
        dot(displacement, cross(b, c)) / determinant,
        dot(displacement, cross(c, a)) / determinant,
        dot(displacement, cross(a, b)) / determinant,
    };
    const std::array<long long, 3> nearest{
        std::llround(fractional[0]), std::llround(fractional[1]), std::llround(fractional[2])};

    Vec3 best = displacement;
    double best_squared = std::numeric_limits<double>::infinity();
    // Check neighboring lattice images around the fractional-coordinate
    // rounding result. This handles skewed conventional simulation cells,
    // where independently wrapping each fractional coordinate is insufficient.
    for (long long i = nearest[0] - 1; i <= nearest[0] + 1; ++i) {
        for (long long j = nearest[1] - 1; j <= nearest[1] + 1; ++j) {
            for (long long k = nearest[2] - 1; k <= nearest[2] + 1; ++k) {
                Vec3 candidate{};
                for (std::size_t dimension = 0; dimension < 3; ++dimension) {
                    candidate[dimension] = displacement[dimension] - static_cast<double>(i) * a[dimension] -
                                           static_cast<double>(j) * b[dimension] - static_cast<double>(k) * c[dimension];
                }
                const double candidate_squared = squared_norm(candidate);
                if (candidate_squared < best_squared) {
                    best = candidate;
                    best_squared = candidate_squared;
                }
            }
        }
    }
    return best;
}

double minimum_image_distance_cell(const Vec3& a, const Vec3& b, const PeriodicCell& cell) {
    const Vec3 displacement{a[0] - b[0], a[1] - b[1], a[2] - b[2]};
    return std::sqrt(squared_norm(minimum_image_displacement(displacement, cell)));
}

double minimum_distance(const Frame& frame, const std::vector<std::size_t>& a, const std::vector<std::size_t>& b) {
    validate_indices(frame, a); validate_indices(frame, b);
    const auto cell = frame_cell(frame);
    double result = std::numeric_limits<double>::infinity();
    for (const auto i : a) for (const auto j : b) if (i != j) result = std::min(result, frame_distance(cell, frame.coordinates[i], frame.coordinates[j]));
    if (!std::isfinite(result)) throw std::invalid_argument("Selections contain no distinct atom pair");
    return result;
}

std::vector<Contact> contacts(const Frame& frame, const std::vector<std::size_t>& a,
                              const std::vector<std::size_t>& b, double cutoff_nm) {
    validate_indices(frame, a); validate_indices(frame, b);
    const auto cell = frame_cell(frame);
    if (cutoff_nm < 0.0) throw std::invalid_argument("Contact cutoff must be non-negative");
    std::vector<Contact> result;
    for (const auto i : a) for (const auto j : b) {
        if (i == j) continue;
        const double d = frame_distance(cell, frame.coordinates[i], frame.coordinates[j]);
        const double tolerance = 1e-12 * std::max(1.0, std::abs(cutoff_nm));
        if (d <= cutoff_nm + tolerance) result.push_back({i, j, d}); // inclusive, robust to unit conversion roundoff
    }
    return result;
}

double radius_of_gyration(const Frame& frame, const std::vector<std::size_t>& selection) {
    return radius_of_gyration(frame, selection, {});
}

double radius_of_gyration(const Frame& frame, const std::vector<std::size_t>& selection,
                          const std::vector<Bond>& bonds) {
    validate_indices(frame, selection);
    std::vector<Vec3> coordinates = frame.coordinates;
    const auto cell = frame_cell(frame);
    if (cell && selection.size() > 1) {
        std::vector<std::vector<std::size_t>> adjacency(frame.coordinates.size());
        for (const auto& bond : bonds) {
            if (bond[0] >= frame.coordinates.size() || bond[1] >= frame.coordinates.size()) {
                throw std::out_of_range("Topology bond index exceeds frame coordinates");
            }
            adjacency[bond[0]].push_back(bond[1]);
            adjacency[bond[1]].push_back(bond[0]);
        }

        std::vector<bool> visited(frame.coordinates.size(), false);
        std::queue<std::size_t> pending;
        const std::size_t root = selection.front();
        visited[root] = true;
        pending.push(root);
        while (!pending.empty()) {
            const std::size_t current = pending.front();
            pending.pop();
            for (const std::size_t neighbor : adjacency[current]) {
                if (visited[neighbor]) continue;
                const Vec3 wrapped_delta{
                    frame.coordinates[neighbor][0] - frame.coordinates[current][0],
                    frame.coordinates[neighbor][1] - frame.coordinates[current][1],
                    frame.coordinates[neighbor][2] - frame.coordinates[current][2],
                };
                const Vec3 delta = minimum_image_displacement(wrapped_delta, *cell);
                for (std::size_t dimension = 0; dimension < 3; ++dimension) {
                    coordinates[neighbor][dimension] = coordinates[current][dimension] + delta[dimension];
                }
                visited[neighbor] = true;
                pending.push(neighbor);
            }
        }
        for (const std::size_t index : selection) {
            if (!visited[index]) {
                throw std::runtime_error(
                    "Periodic RG selection spans atoms that are not connected by topology bonds; "
                    "provide PDB CONECT records for molecule unwrapping");
            }
        }
    }
    Vec3 center{};
    for (const auto i : selection) for (std::size_t d = 0; d < 3; ++d) center[d] += coordinates[i][d];
    for (double& value : center) value /= static_cast<double>(selection.size());
    double squared = 0.0;
    for (const auto i : selection) {
        const double d = distance(coordinates[i], center);
        squared += d * d;
    }
    return std::sqrt(squared / static_cast<double>(selection.size()));
}

} // namespace ensembleql
