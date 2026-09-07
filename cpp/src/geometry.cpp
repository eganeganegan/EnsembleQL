#include "ensembleql/geometry.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <functional>
#include <limits>
#include <optional>
#include <queue>
#include <set>
#include <stdexcept>
#include <tuple>
#include <unordered_map>
#include <utility>

namespace ensembleql {
namespace {
constexpr std::size_t spatial_hash_pair_threshold = 4096;

void validate_indices(const Frame& frame, const std::vector<std::size_t>& indices) {
    if (indices.empty()) throw std::invalid_argument("Geometry selection is empty");
    for (const auto index : indices) {
        if (index >= frame.coordinates.size()) throw std::out_of_range("Atom index exceeds frame coordinates");
        for (const double coordinate : frame.coordinates[index]) {
            if (!std::isfinite(coordinate)) {
                throw std::invalid_argument("Selected atom coordinates must be finite");
            }
        }
    }
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

bool same_cell(const std::optional<PeriodicCell>& first,
               const std::optional<PeriodicCell>& second) {
    if (first.has_value() != second.has_value()) return false;
    if (!first) return true;
    for (std::size_t vector = 0; vector < 3; ++vector) {
        for (std::size_t dimension = 0; dimension < 3; ++dimension) {
            if (first->vectors_nm[vector][dimension] !=
                second->vectors_nm[vector][dimension]) return false;
        }
    }
    return true;
}

std::optional<Vec3> axis_aligned_box(const Frame& frame) {
    if (!frame.cell_nm) {
        if (frame.box_nm) validate_box(*frame.box_nm);
        return frame.box_nm;
    }
    validate_cell(*frame.cell_nm);
    constexpr double zero_tolerance = 1e-14;
    Vec3 box{};
    for (std::size_t vector = 0; vector < 3; ++vector) {
        for (std::size_t dimension = 0; dimension < 3; ++dimension) {
            const double value = frame.cell_nm->vectors_nm[vector][dimension];
            if (vector == dimension) {
                if (value <= 0.0) return std::nullopt;
                box[dimension] = value;
            } else if (std::abs(value) > zero_tolerance) {
                return std::nullopt;
            }
        }
    }
    return box;
}

using CellKey = std::array<long long, 3>;

struct CellKeyHash {
    std::size_t operator()(const CellKey& key) const noexcept {
        std::size_t result = 0;
        for (const long long value : key) {
            const std::size_t hash = std::hash<long long>{}(value);
            result ^= hash + static_cast<std::size_t>(0x9e3779b9U) + (result << 6U) + (result >> 2U);
        }
        return result;
    }
};

struct SpatialGrid {
    double cell_width_nm{};
    std::optional<Vec3> box_nm;
    std::optional<PeriodicCell> cell_nm;
    std::array<long long, 3> periodic_bins{};
    std::unordered_map<CellKey, std::vector<std::size_t>, CellKeyHash> cells;
};

bool pair_count_reaches_threshold(std::size_t a_size, std::size_t b_size) {
    if (b_size == 0) return false;
    if (b_size >= spatial_hash_pair_threshold) return a_size != 0;
    return a_size >= (spatial_hash_pair_threshold + b_size - 1) / b_size;
}

std::optional<CellKey> nonperiodic_key(const Vec3& coordinate, double width) {
    CellKey key{};
    constexpr double key_margin = 4096.0;
    constexpr double minimum_key = static_cast<double>(std::numeric_limits<long long>::lowest()) + key_margin;
    constexpr double maximum_key = static_cast<double>(std::numeric_limits<long long>::max()) - key_margin;
    for (std::size_t dimension = 0; dimension < 3; ++dimension) {
        const double scaled = std::floor(coordinate[dimension] / width);
        if (!std::isfinite(scaled) || scaled < minimum_key || scaled > maximum_key) return std::nullopt;
        key[dimension] = static_cast<long long>(scaled);
    }
    return key;
}

Vec3 fractional_coordinate(const Vec3& coordinate, const PeriodicCell& cell) {
    const auto& a = cell.vectors_nm[0];
    const auto& b = cell.vectors_nm[1];
    const auto& c = cell.vectors_nm[2];
    const double determinant = dot(a, cross(b, c));
    return {dot(coordinate, cross(b, c)) / determinant,
            dot(coordinate, cross(c, a)) / determinant,
            dot(coordinate, cross(a, b)) / determinant};
}

CellKey periodic_key(const Vec3& coordinate, const SpatialGrid& grid) {
    CellKey key{};
    if (grid.cell_nm) {
        const Vec3 fractional = fractional_coordinate(coordinate, *grid.cell_nm);
        for (std::size_t dimension = 0; dimension < 3; ++dimension) {
            double wrapped = fractional[dimension] - std::floor(fractional[dimension]);
            if (wrapped >= 1.0) wrapped = 0.0;
            key[dimension] = std::min(
                grid.periodic_bins[dimension] - 1,
                static_cast<long long>(std::floor(
                    wrapped * static_cast<double>(grid.periodic_bins[dimension]))));
        }
        return key;
    }
    for (std::size_t dimension = 0; dimension < 3; ++dimension) {
        const double length = (*grid.box_nm)[dimension];
        double wrapped = std::fmod(coordinate[dimension], length);
        if (wrapped < 0.0) wrapped += length;
        const double bin_width = length / static_cast<double>(grid.periodic_bins[dimension]);
        key[dimension] = std::min(grid.periodic_bins[dimension] - 1,
                                  static_cast<long long>(std::floor(wrapped / bin_width)));
    }
    return key;
}

std::optional<SpatialGrid> make_spatial_grid(const Frame& frame,
                                             const std::vector<std::size_t>& selection,
                                             double cell_width_nm) {
    SpatialGrid grid;
    grid.cell_width_nm = cell_width_nm;
    grid.box_nm = axis_aligned_box(frame);
    if (frame.cell_nm && !grid.box_nm) grid.cell_nm = frame.cell_nm;

    if (grid.box_nm || grid.cell_nm) {
        for (std::size_t dimension = 0; dimension < 3; ++dimension) {
            double ratio{};
            if (grid.box_nm) {
                ratio = std::floor((*grid.box_nm)[dimension] / cell_width_nm);
            } else {
                const auto& cell = *grid.cell_nm;
                const Vec3 reciprocal = dimension == 0 ? cross(cell.vectors_nm[1], cell.vectors_nm[2]) :
                    dimension == 1 ? cross(cell.vectors_nm[2], cell.vectors_nm[0]) :
                                     cross(cell.vectors_nm[0], cell.vectors_nm[1]);
                const double determinant = std::abs(dot(
                    cell.vectors_nm[0], cross(cell.vectors_nm[1], cell.vectors_nm[2])));
                const double fractional_reach = cell_width_nm * std::sqrt(squared_norm(reciprocal)) /
                                                determinant;
                ratio = std::floor(1.0 / fractional_reach);
            }
            constexpr double maximum_bins =
                static_cast<double>(std::numeric_limits<long long>::max() / 2);
            if (!std::isfinite(ratio) || ratio > maximum_bins) {
                return std::nullopt;
            }
            grid.periodic_bins[dimension] = std::max(1LL, static_cast<long long>(ratio));
        }
    } else {
        for (const std::size_t index : selection) {
            if (!nonperiodic_key(frame.coordinates[index], cell_width_nm)) return std::nullopt;
        }
    }

    grid.cells.reserve(selection.size());
    for (std::size_t position = 0; position < selection.size(); ++position) {
        const Vec3& coordinate = frame.coordinates[selection[position]];
        const CellKey key = (grid.box_nm || grid.cell_nm)
                                ? periodic_key(coordinate, grid)
                                : *nonperiodic_key(coordinate, cell_width_nm);
        grid.cells[key].push_back(position);
    }
    return grid;
}

long long wrapped_bin(long long value, long long count) {
    value %= count;
    return value < 0 ? value + count : value;
}

std::vector<std::size_t> spatial_candidates(const SpatialGrid& grid, const Vec3& coordinate) {
    const auto center = (grid.box_nm || grid.cell_nm)
                            ? std::optional<CellKey>{periodic_key(coordinate, grid)}
                            : nonperiodic_key(coordinate, grid.cell_width_nm);
    if (!center) return {};

    std::array<CellKey, 27> neighbor_keys{};
    std::size_t neighbor_count = 0;
    for (long long x = -1; x <= 1; ++x) {
        for (long long y = -1; y <= 1; ++y) {
            for (long long z = -1; z <= 1; ++z) {
                CellKey key{(*center)[0] + x, (*center)[1] + y, (*center)[2] + z};
                if (grid.box_nm || grid.cell_nm) {
                    for (std::size_t dimension = 0; dimension < 3; ++dimension) {
                        key[dimension] = wrapped_bin(key[dimension], grid.periodic_bins[dimension]);
                    }
                }
                neighbor_keys[neighbor_count++] = key;
            }
        }
    }
    std::sort(neighbor_keys.begin(), neighbor_keys.end());
    const auto unique_end = std::unique(neighbor_keys.begin(), neighbor_keys.end());

    std::vector<std::size_t> candidates;
    for (auto key = neighbor_keys.begin(); key != unique_end; ++key) {
        const auto found = grid.cells.find(*key);
        if (found != grid.cells.end()) {
            candidates.insert(candidates.end(), found->second.begin(), found->second.end());
        }
    }
    std::sort(candidates.begin(), candidates.end());
    return candidates;
}

template <typename Visitor>
bool visit_contacts_brute_force(const Frame& frame,
                                const std::vector<std::size_t>& a,
                                const std::vector<std::size_t>& b,
                                const std::optional<PeriodicCell>& cell,
                                double threshold_nm,
                                Visitor&& visitor) {
    for (const std::size_t i : a) {
        for (const std::size_t j : b) {
            if (i == j) continue;
            const double separation = frame_distance(cell, frame.coordinates[i], frame.coordinates[j]);
            if (separation <= threshold_nm && visitor(i, j, separation)) return true;
        }
    }
    return false;
}

template <typename Visitor>
bool visit_contacts(const Frame& frame,
                    const std::vector<std::size_t>& a,
                    const std::vector<std::size_t>& b,
                    double cutoff_nm,
                    Visitor&& visitor) {
    if (!std::isfinite(cutoff_nm) || cutoff_nm < 0.0) {
        throw std::invalid_argument("Contact cutoff must be finite and non-negative");
    }
    const auto cell = frame_cell(frame);
    const double tolerance = 1e-12 * std::max(1.0, std::abs(cutoff_nm));
    const double threshold_nm = cutoff_nm > std::numeric_limits<double>::max() - tolerance
                                    ? std::numeric_limits<double>::max()
                                    : cutoff_nm + tolerance;
    if (threshold_nm <= 0.0 || !pair_count_reaches_threshold(a.size(), b.size())) {
        return visit_contacts_brute_force(frame, a, b, cell, threshold_nm,
                                          std::forward<Visitor>(visitor));
    }

    auto grid = make_spatial_grid(frame, b, threshold_nm);
    if (!grid) {
        return visit_contacts_brute_force(frame, a, b, cell, threshold_nm,
                                          std::forward<Visitor>(visitor));
    }
    if (!grid->box_nm && !grid->cell_nm) {
        for (const std::size_t index : a) {
            if (!nonperiodic_key(frame.coordinates[index], threshold_nm)) {
                return visit_contacts_brute_force(frame, a, b, cell, threshold_nm,
                                                  std::forward<Visitor>(visitor));
            }
        }
    }
    for (const std::size_t i : a) {
        const auto candidates = spatial_candidates(*grid, frame.coordinates[i]);
        for (const std::size_t position : candidates) {
            const std::size_t j = b[position];
            if (i == j) continue;
            const double separation = frame_distance(cell, frame.coordinates[i], frame.coordinates[j]);
            if (separation <= threshold_nm && visitor(i, j, separation)) return true;
        }
    }
    return false;
}

bool has_distinct_pair(const std::vector<std::size_t>& a, const std::vector<std::size_t>& b) {
    for (const std::size_t j : b) if (a.front() != j) return true;
    for (const std::size_t i : a) if (i != a.front()) return true;
    return false;
}

std::vector<Vec3> coordinates_for_rg(const Frame& frame, const std::vector<std::size_t>& selection,
                                     const std::vector<Bond>& bonds) {
    std::vector<Vec3> coordinates = frame.coordinates;
    const auto cell = frame_cell(frame);
    if (!cell || selection.size() <= 1) return coordinates;

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
    return coordinates;
}

double rg_from_coordinates(const std::vector<Vec3>& coordinates,
                           const std::vector<std::size_t>& selection,
                           const std::vector<double>* masses_da) {
    Vec3 center{};
    double total_weight = 0.0;
    for (const auto index : selection) {
        const double weight = masses_da ? (*masses_da)[index] : 1.0;
        if (!std::isfinite(weight) || weight <= 0.0) {
            throw std::runtime_error("Mass-weighted RG requires a known positive mass for every selected atom");
        }
        total_weight += weight;
        for (std::size_t dimension = 0; dimension < 3; ++dimension) {
            center[dimension] += weight * coordinates[index][dimension];
        }
    }
    for (double& value : center) value /= total_weight;
    double squared = 0.0;
    for (const auto index : selection) {
        const double weight = masses_da ? (*masses_da)[index] : 1.0;
        const double displacement = distance(coordinates[index], center);
        squared += weight * displacement * displacement;
    }
    return std::sqrt(squared / total_weight);
}
} // namespace

NeighborList::NeighborList(std::vector<std::size_t> selection_a,
                           std::vector<std::size_t> selection_b,
                           double cutoff_nm, double skin_nm)
    : selection_a_(std::move(selection_a)), selection_b_(std::move(selection_b)),
      cutoff_nm_(cutoff_nm), skin_nm_(skin_nm) {
    if (!std::isfinite(cutoff_nm_) || cutoff_nm_ < 0.0) {
        throw std::invalid_argument("Neighbor-list cutoff must be finite and non-negative");
    }
    if (!std::isfinite(skin_nm_) || skin_nm_ <= 0.0) {
        throw std::invalid_argument("Neighbor-list skin must be finite and positive");
    }
    if (cutoff_nm_ > std::numeric_limits<double>::max() - skin_nm_) {
        throw std::invalid_argument("Neighbor-list cutoff plus skin must be finite");
    }
    tracked_atoms_ = selection_a_;
    tracked_atoms_.insert(tracked_atoms_.end(), selection_b_.begin(), selection_b_.end());
    std::sort(tracked_atoms_.begin(), tracked_atoms_.end());
    tracked_atoms_.erase(std::unique(tracked_atoms_.begin(), tracked_atoms_.end()),
                         tracked_atoms_.end());
}

bool NeighborList::needs_rebuild(const Frame& frame) const {
    if (!initialized_ || frame.coordinates.size() != reference_coordinate_count_) return true;
    const auto cell = frame_cell(frame);
    if (!same_cell(cell, reference_cell_)) return true;
    const double displacement_limit = skin_nm_ * 0.5;
    for (std::size_t position = 0; position < tracked_atoms_.size(); ++position) {
        const std::size_t index = tracked_atoms_[position];
        if (index >= frame.coordinates.size()) return true;
        for (const double coordinate : frame.coordinates[index]) {
            if (!std::isfinite(coordinate)) {
                throw std::invalid_argument("Selected atom coordinates must be finite");
            }
        }
        if (frame_distance(cell, frame.coordinates[index], reference_coordinates_[position]) >
            displacement_limit) return true;
    }
    return false;
}

void NeighborList::rebuild(const Frame& frame) {
    validate_indices(frame, selection_a_);
    validate_indices(frame, selection_b_);
    candidate_pairs_.clear();
    const auto candidates = ensembleql::contacts(
        frame, selection_a_, selection_b_, cutoff_nm_ + skin_nm_);
    candidate_pairs_.reserve(candidates.size());
    for (const Contact& candidate : candidates) {
        candidate_pairs_.emplace_back(candidate.atom_a, candidate.atom_b);
    }
    reference_coordinates_.clear();
    reference_coordinates_.reserve(tracked_atoms_.size());
    for (const std::size_t index : tracked_atoms_) {
        reference_coordinates_.push_back(frame.coordinates[index]);
    }
    reference_cell_ = frame_cell(frame);
    reference_coordinate_count_ = frame.coordinates.size();
    initialized_ = true;
    ++rebuild_count_;
}

void NeighborList::update(const Frame& frame) {
    if (needs_rebuild(frame)) rebuild(frame);
}

std::vector<Contact> NeighborList::contacts(const Frame& frame) {
    update(frame);
    const auto cell = frame_cell(frame);
    const double tolerance = 1e-12 * std::max(1.0, std::abs(cutoff_nm_));
    const double threshold = cutoff_nm_ > std::numeric_limits<double>::max() - tolerance
                                 ? std::numeric_limits<double>::max()
                                 : cutoff_nm_ + tolerance;
    std::vector<Contact> result;
    if (!cell) {
        std::vector<double> separations(candidate_pairs_.size());
#ifdef ENSEMBLEQL_HAS_OPENMP
#pragma omp simd
#endif
        for (std::ptrdiff_t position = 0;
             position < static_cast<std::ptrdiff_t>(candidate_pairs_.size()); ++position) {
            const auto& pair = candidate_pairs_[static_cast<std::size_t>(position)];
            const Vec3& first = frame.coordinates[pair.first];
            const Vec3& second = frame.coordinates[pair.second];
            const double dx = first[0] - second[0];
            const double dy = first[1] - second[1];
            const double dz = first[2] - second[2];
            separations[static_cast<std::size_t>(position)] =
                std::sqrt(dx * dx + dy * dy + dz * dz);
        }
        for (std::size_t position = 0; position < candidate_pairs_.size(); ++position) {
            if (separations[position] <= threshold) {
                result.push_back({candidate_pairs_[position].first,
                                  candidate_pairs_[position].second,
                                  separations[position]});
            }
        }
        return result;
    }
    for (const auto& [first, second] : candidate_pairs_) {
        const double separation = frame_distance(
            cell, frame.coordinates[first], frame.coordinates[second]);
        if (separation <= threshold) result.push_back({first, second, separation});
    }
    return result;
}

bool NeighborList::has_contact(const Frame& frame) {
    update(frame);
    bool distinct_pair = false;
    for (const std::size_t first : selection_a_) {
        if (std::any_of(selection_b_.begin(), selection_b_.end(),
                        [first](std::size_t second) { return first != second; })) {
            distinct_pair = true;
            break;
        }
    }
    if (!distinct_pair) {
        throw std::invalid_argument("Selections contain no distinct atom pair");
    }
    const auto cell = frame_cell(frame);
    const double tolerance = 1e-12 * std::max(1.0, std::abs(cutoff_nm_));
    const double threshold = cutoff_nm_ > std::numeric_limits<double>::max() - tolerance
                                 ? std::numeric_limits<double>::max()
                                 : cutoff_nm_ + tolerance;
    for (const auto& [first, second] : candidate_pairs_) {
        if (frame_distance(cell, frame.coordinates[first], frame.coordinates[second]) <=
            threshold) return true;
    }
    return false;
}

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
    std::vector<Contact> result;
    visit_contacts(frame, a, b, cutoff_nm,
                   [&result](std::size_t i, std::size_t j, double separation) {
                       result.push_back({i, j, separation});
                       return false;
                   });
    return result;
}

bool has_contact(const Frame& frame, const std::vector<std::size_t>& a,
                 const std::vector<std::size_t>& b, double cutoff_nm) {
    validate_indices(frame, a); validate_indices(frame, b);
    if (!has_distinct_pair(a, b)) {
        throw std::invalid_argument("Selections contain no distinct atom pair");
    }
    return visit_contacts(frame, a, b, cutoff_nm,
                          [](std::size_t, std::size_t, double) { return true; });
}

std::size_t contact_count(const Frame& frame, const std::vector<std::size_t>& a,
                          const std::vector<std::size_t>& b, double cutoff_nm,
                          ContactMode mode, const Topology& topology) {
    const auto atom_contacts = contacts(frame, a, b, cutoff_nm);
    std::set<Bond> unique_atom_contacts;
    for (const Contact& contact : atom_contacts) {
        Bond pair{contact.atom_a, contact.atom_b};
        if (pair[1] < pair[0]) std::swap(pair[0], pair[1]);
        unique_atom_contacts.insert(pair);
    }
    if (mode == ContactMode::Atom) return unique_atom_contacts.size();
    if (topology.size() != frame.coordinates.size()) {
        throw std::invalid_argument("Residue contact counting requires topology and frame atom counts to match");
    }
    using Residue = std::pair<char, int>;
    using ResiduePair = std::pair<Residue, Residue>;
    std::set<ResiduePair> residue_contacts;
    for (const Bond& contact : unique_atom_contacts) {
        Residue first{topology.atoms()[contact[0]].chain,
                      topology.atoms()[contact[0]].residue_index};
        Residue second{topology.atoms()[contact[1]].chain,
                       topology.atoms()[contact[1]].residue_index};
        if (first == second) continue;
        if (second < first) std::swap(first, second);
        residue_contacts.emplace(first, second);
    }
    return residue_contacts.size();
}

double radius_of_gyration(const Frame& frame, const std::vector<std::size_t>& selection) {
    return radius_of_gyration(frame, selection, {});
}

double radius_of_gyration(const Frame& frame, const std::vector<std::size_t>& selection,
                          const std::vector<Bond>& bonds) {
    validate_indices(frame, selection);
    return rg_from_coordinates(coordinates_for_rg(frame, selection, bonds), selection, nullptr);
}

double mass_weighted_radius_of_gyration(const Frame& frame,
                                        const std::vector<std::size_t>& selection,
                                        const std::vector<Bond>& bonds,
                                        const std::vector<double>& masses_da) {
    validate_indices(frame, selection);
    if (masses_da.size() != frame.coordinates.size()) {
        throw std::invalid_argument("Mass-weighted RG requires one atom mass per frame coordinate");
    }
    return rg_from_coordinates(coordinates_for_rg(frame, selection, bonds), selection, &masses_da);
}

} // namespace ensembleql
