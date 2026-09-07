#include "ensembleql/molecular.hpp"

#include "ensembleql/geometry.hpp"

#include <algorithm>
#include <array>
#include <cctype>
#include <cmath>
#include <limits>
#include <map>
#include <optional>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>

namespace ensembleql {
namespace {
constexpr double pi = 3.14159265358979323846;

void validate_selection(const Frame& frame, const std::vector<std::size_t>& selection,
                        const char* name) {
    if (selection.empty()) throw std::invalid_argument(std::string(name) + " selection is empty");
    for (const std::size_t index : selection) {
        if (index >= frame.coordinates.size()) {
            throw std::out_of_range(std::string(name) + " atom index exceeds frame coordinates");
        }
        for (const double value : frame.coordinates[index]) {
            if (!std::isfinite(value)) throw std::invalid_argument("Molecular coordinates must be finite");
        }
    }
}

PeriodicCell orthorhombic_cell(const Vec3& box) {
    for (const double length : box) {
        if (!std::isfinite(length) || length <= 0.0) {
            throw std::invalid_argument("Periodic box lengths must be finite and positive");
        }
    }
    PeriodicCell cell;
    cell.vectors_nm = {Vec3{box[0], 0.0, 0.0}, Vec3{0.0, box[1], 0.0},
                       Vec3{0.0, 0.0, box[2]}};
    return cell;
}

std::optional<PeriodicCell> frame_cell(const Frame& frame) {
    if (frame.cell_nm) return frame.cell_nm;
    if (frame.box_nm) return orthorhombic_cell(*frame.box_nm);
    return std::nullopt;
}

Vec3 subtract(const Vec3& a, const Vec3& b) {
    return {a[0] - b[0], a[1] - b[1], a[2] - b[2]};
}

Vec3 add(const Vec3& a, const Vec3& b) {
    return {a[0] + b[0], a[1] + b[1], a[2] + b[2]};
}

Vec3 scale(const Vec3& value, double factor) {
    return {value[0] * factor, value[1] * factor, value[2] * factor};
}

double dot(const Vec3& a, const Vec3& b) {
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

Vec3 cross(const Vec3& a, const Vec3& b) {
    return {a[1] * b[2] - a[2] * b[1],
            a[2] * b[0] - a[0] * b[2],
            a[0] * b[1] - a[1] * b[0]};
}

double norm(const Vec3& value) {
    return std::sqrt(dot(value, value));
}

Vec3 normalized(const Vec3& value) {
    const double length = norm(value);
    if (!std::isfinite(length) || length <= 1e-15) {
        throw std::invalid_argument("Molecular geometry contains a zero-length vector");
    }
    return scale(value, 1.0 / length);
}

Vec3 displacement(const Frame& frame, const Vec3& from, const Vec3& to) {
    const Vec3 raw = subtract(to, from);
    const auto cell = frame_cell(frame);
    return cell ? minimum_image_displacement(raw, *cell) : raw;
}

Vec3 displacement(const Frame& frame, std::size_t from, std::size_t to) {
    return displacement(frame, frame.coordinates[from], frame.coordinates[to]);
}

std::vector<Vec3> unwrapped(const Frame& frame, const std::vector<std::size_t>& selection) {
    validate_selection(frame, selection, "Geometry");
    std::vector<Vec3> result;
    result.reserve(selection.size());
    const Vec3 anchor = frame.coordinates[selection.front()];
    result.push_back(anchor);
    for (std::size_t position = 1; position < selection.size(); ++position) {
        result.push_back(add(anchor, displacement(frame, anchor, frame.coordinates[selection[position]])));
    }
    return result;
}

Vec3 centroid(const std::vector<Vec3>& coordinates) {
    if (coordinates.empty()) throw std::invalid_argument("Centroid requires coordinates");
    Vec3 result{};
    for (const auto& coordinate : coordinates) result = add(result, coordinate);
    return scale(result, 1.0 / static_cast<double>(coordinates.size()));
}

double degrees(double radians) {
    return radians * 180.0 / pi;
}

double vector_angle(const Vec3& first, const Vec3& second, bool unoriented = false) {
    double cosine = dot(normalized(first), normalized(second));
    if (unoriented) cosine = std::abs(cosine);
    cosine = std::clamp(cosine, -1.0, 1.0);
    return degrees(std::acos(cosine));
}

std::string uppercase(std::string value) {
    std::transform(value.begin(), value.end(), value.begin(), [](unsigned char character) {
        return static_cast<char>(std::toupper(character));
    });
    return value;
}

const std::map<std::string, double>& vdw_radii_nm() {
    static const std::map<std::string, double> radii{
        {"H", 0.120}, {"C", 0.170}, {"N", 0.155}, {"O", 0.152},
        {"F", 0.147}, {"SI", 0.210}, {"P", 0.180}, {"S", 0.180}, {"CL", 0.175},
        {"BR", 0.185}, {"I", 0.198},
    };
    return radii;
}

double vdw_radius_nm(const std::string& element) {
    const auto& radii = vdw_radii_nm();
    const auto found = radii.find(uppercase(element));
    if (found == radii.end()) {
        throw std::runtime_error("SASA requires a van der Waals radius for element '" + element + "'");
    }
    return found->second;
}

std::vector<Vec3> fibonacci_sphere(std::size_t count) {
    std::vector<Vec3> points;
    points.reserve(count);
    constexpr double golden_angle = 2.39996322972865332;
    for (std::size_t index = 0; index < count; ++index) {
        const double y = 1.0 - 2.0 * (static_cast<double>(index) + 0.5) /
                                   static_cast<double>(count);
        const double radius = std::sqrt(std::max(0.0, 1.0 - y * y));
        const double phi = golden_angle * static_cast<double>(index);
        points.push_back({std::cos(phi) * radius, y, std::sin(phi) * radius});
    }
    return points;
}

Vec3 ring_normal(const std::vector<Vec3>& coordinates) {
    const Vec3 center = centroid(coordinates);
    Vec3 normal{};
    for (std::size_t index = 0; index < coordinates.size(); ++index) {
        const Vec3 first = subtract(coordinates[index], center);
        const Vec3 second = subtract(coordinates[(index + 1) % coordinates.size()], center);
        normal = add(normal, cross(first, second));
    }
    return normalized(normal);
}

using ResidueKey = std::pair<char, int>;
struct BackboneAtoms {
    std::optional<std::size_t> n;
    std::optional<std::size_t> ca;
    std::optional<std::size_t> c;
};

std::array<double, 4> largest_eigenvector(std::array<std::array<double, 4>, 4> matrix) {
    std::array<std::array<double, 4>, 4> vectors{};
    for (std::size_t i = 0; i < 4; ++i) vectors[i][i] = 1.0;
    for (std::size_t iteration = 0; iteration < 64; ++iteration) {
        std::size_t p = 0, q = 1;
        double largest = std::abs(matrix[p][q]);
        for (std::size_t i = 0; i < 4; ++i) {
            for (std::size_t j = i + 1; j < 4; ++j) {
                if (std::abs(matrix[i][j]) > largest) {
                    largest = std::abs(matrix[i][j]); p = i; q = j;
                }
            }
        }
        if (largest < 1e-14) break;
        const double angle = 0.5 * std::atan2(2.0 * matrix[p][q], matrix[q][q] - matrix[p][p]);
        const double cosine = std::cos(angle), sine = std::sin(angle);
        for (std::size_t k = 0; k < 4; ++k) {
            const double mkp = matrix[k][p], mkq = matrix[k][q];
            matrix[k][p] = cosine * mkp - sine * mkq;
            matrix[k][q] = sine * mkp + cosine * mkq;
        }
        for (std::size_t k = 0; k < 4; ++k) {
            const double mpk = matrix[p][k], mqk = matrix[q][k];
            matrix[p][k] = cosine * mpk - sine * mqk;
            matrix[q][k] = sine * mpk + cosine * mqk;
        }
        for (std::size_t k = 0; k < 4; ++k) {
            const double vkp = vectors[k][p], vkq = vectors[k][q];
            vectors[k][p] = cosine * vkp - sine * vkq;
            vectors[k][q] = sine * vkp + cosine * vkq;
        }
    }
    std::size_t largest = 0;
    for (std::size_t i = 1; i < 4; ++i) if (matrix[i][i] > matrix[largest][largest]) largest = i;
    std::array<double, 4> result{};
    for (std::size_t i = 0; i < 4; ++i) result[i] = vectors[i][largest];
    return result;
}

Vec3 rotate(const Vec3& point, const std::array<double, 4>& quaternion) {
    const double w = quaternion[0], x = quaternion[1], y = quaternion[2], z = quaternion[3];
    return {
        (1 - 2 * (y * y + z * z)) * point[0] + 2 * (x * y - z * w) * point[1] + 2 * (x * z + y * w) * point[2],
        2 * (x * y + z * w) * point[0] + (1 - 2 * (x * x + z * z)) * point[1] + 2 * (y * z - x * w) * point[2],
        2 * (x * z - y * w) * point[0] + 2 * (y * z + x * w) * point[1] + (1 - 2 * (x * x + y * y)) * point[2],
    };
}
} // namespace

bool sasa_element_supported(const std::string& element) {
    return vdw_radii_nm().contains(uppercase(element));
}

double angle_degrees(const Frame& frame, std::size_t first,
                     std::size_t vertex, std::size_t third) {
    validate_selection(frame, {first, vertex, third}, "Angle");
    return vector_angle(displacement(frame, vertex, first), displacement(frame, vertex, third));
}

double dihedral_degrees(const Frame& frame, std::size_t first, std::size_t second,
                        std::size_t third, std::size_t fourth) {
    validate_selection(frame, {first, second, third, fourth}, "Dihedral");
    const Vec3 b1 = displacement(frame, first, second);
    const Vec3 b2 = displacement(frame, second, third);
    const Vec3 b3 = displacement(frame, third, fourth);
    const Vec3 n1 = normalized(cross(b1, b2));
    const Vec3 n2 = normalized(cross(b2, b3));
    const Vec3 axis = normalized(b2);
    return degrees(std::atan2(dot(cross(n1, n2), axis), dot(n1, n2)));
}

bool hydrogen_bond(const Frame& frame, const std::vector<std::size_t>& donors,
                   const std::vector<std::size_t>& acceptors, const Topology& topology,
                   double distance_cutoff_nm, double minimum_angle_degrees) {
    validate_selection(frame, donors, "Hydrogen-bond donor");
    validate_selection(frame, acceptors, "Hydrogen-bond acceptor");
    if (topology.size() != frame.coordinates.size()) throw std::invalid_argument("Hydrogen bonds require matching topology and frame atom counts");
    if (!std::isfinite(distance_cutoff_nm) || distance_cutoff_nm < 0.0 ||
        !std::isfinite(minimum_angle_degrees) || minimum_angle_degrees < 0.0 ||
        minimum_angle_degrees > 180.0) {
        throw std::invalid_argument("Invalid hydrogen-bond distance or angle criterion");
    }
    for (const std::size_t donor : donors) {
        std::vector<std::size_t> hydrogens;
        for (const Bond& bond : topology.bonds()) {
            const std::size_t neighbor = bond[0] == donor ? bond[1] : bond[1] == donor ? bond[0] : topology.size();
            if (neighbor < topology.size() && uppercase(topology.atoms()[neighbor].element) == "H") {
                hydrogens.push_back(neighbor);
            }
        }
        for (const std::size_t acceptor : acceptors) {
            if (donor == acceptor || norm(displacement(frame, donor, acceptor)) >
                                         distance_cutoff_nm + 1e-12) continue;
            for (const std::size_t hydrogen : hydrogens) {
                if (angle_degrees(frame, donor, hydrogen, acceptor) + 1e-12 >= minimum_angle_degrees) return true;
            }
        }
    }
    return false;
}

double helix_fraction(const Frame& frame, const std::vector<std::size_t>& selection,
                      const Topology& topology) {
    validate_selection(frame, selection, "HELIX");
    if (topology.size() != frame.coordinates.size()) throw std::invalid_argument("HELIX requires matching topology and frame atom counts");
    std::map<ResidueKey, BackboneAtoms> residues;
    for (const std::size_t index : selection) {
        const Atom& atom = topology.atoms()[index];
        auto& backbone = residues[{atom.chain, atom.residue_index}];
        const std::string name = uppercase(atom.name);
        if (name == "N") backbone.n = index;
        else if (name == "CA") backbone.ca = index;
        else if (name == "C") backbone.c = index;
    }
    std::vector<std::pair<ResidueKey, BackboneAtoms>> ordered(residues.begin(), residues.end());
    std::size_t evaluable = 0, helical = 0;
    for (std::size_t index = 1; index + 1 < ordered.size(); ++index) {
        if (ordered[index - 1].first.first != ordered[index].first.first ||
            ordered[index].first.first != ordered[index + 1].first.first ||
            ordered[index - 1].first.second + 1 != ordered[index].first.second ||
            ordered[index].first.second + 1 != ordered[index + 1].first.second) continue;
        const auto& previous = ordered[index - 1].second;
        const auto& current = ordered[index].second;
        const auto& next = ordered[index + 1].second;
        if (!previous.c || !current.n || !current.ca || !current.c || !next.n) continue;
        const double phi = dihedral_degrees(frame, *previous.c, *current.n, *current.ca, *current.c);
        const double psi = dihedral_degrees(frame, *current.n, *current.ca, *current.c, *next.n);
        ++evaluable;
        if (phi >= -100.0 && phi <= -30.0 && psi >= -80.0 && psi <= -5.0) ++helical;
    }
    if (evaluable == 0) throw std::runtime_error("HELIX requires at least one residue with complete neighboring N-CA-C backbone atoms");
    return static_cast<double>(helical) / static_cast<double>(evaluable);
}

double solvent_accessible_surface_area(const Frame& frame,
                                       const std::vector<std::size_t>& selection,
                                       const Topology& topology,
                                       double probe_radius_nm, std::size_t sphere_points) {
    validate_selection(frame, selection, "SASA");
    if (topology.size() != frame.coordinates.size()) throw std::invalid_argument("SASA requires matching topology and frame atom counts");
    if (!std::isfinite(probe_radius_nm) || probe_radius_nm < 0.0 || sphere_points < 6) {
        throw std::invalid_argument("SASA requires a non-negative probe radius and at least six sphere points");
    }
    const auto points = fibonacci_sphere(sphere_points);
    std::vector<double> expanded_radii;
    expanded_radii.reserve(topology.size());
    for (const Atom& atom : topology.atoms()) {
        expanded_radii.push_back(vdw_radius_nm(atom.element) + probe_radius_nm);
    }
    double area = 0.0;
#ifdef ENSEMBLEQL_HAS_OPENMP
#pragma omp parallel for reduction(+ : area) schedule(static)
#endif
    for (std::ptrdiff_t position = 0;
         position < static_cast<std::ptrdiff_t>(selection.size()); ++position) {
        const std::size_t atom = selection[static_cast<std::size_t>(position)];
        const double radius = expanded_radii[atom];
        std::size_t accessible = 0;
        for (const Vec3& direction : points) {
            const Vec3 surface = add(frame.coordinates[atom], scale(direction, radius));
            bool occluded = false;
            for (std::size_t other = 0; other < topology.size(); ++other) {
                if (other == atom) continue;
                const double other_radius = expanded_radii[other];
                if (norm(displacement(frame, surface, frame.coordinates[other])) < other_radius) {
                    occluded = true; break;
                }
            }
            if (!occluded) ++accessible;
        }
        area += 4.0 * pi * radius * radius * static_cast<double>(accessible) /
                static_cast<double>(sphere_points);
    }
    return area;
}

double aligned_rmsd(const Frame& frame, const std::vector<std::size_t>& selection,
                    const std::vector<Vec3>& reference, bool align) {
    validate_selection(frame, selection, "RMSD");
    if (reference.size() != selection.size()) throw std::invalid_argument("RMSD reference size does not match selection");
    auto current = unwrapped(frame, selection);
    std::vector<Vec3> target = reference;
    const Vec3 reference_anchor = target.front();
    for (std::size_t atom = 1; atom < target.size(); ++atom) {
        target[atom] = add(reference_anchor, displacement(frame, reference_anchor, target[atom]));
    }
    if (align) {
        const Vec3 current_center = centroid(current), target_center = centroid(target);
        for (Vec3& point : current) point = subtract(point, current_center);
        for (Vec3& point : target) point = subtract(point, target_center);
        std::array<std::array<double, 3>, 3> covariance{};
        for (std::size_t atom = 0; atom < current.size(); ++atom) {
            for (std::size_t i = 0; i < 3; ++i) {
                for (std::size_t j = 0; j < 3; ++j) covariance[i][j] += current[atom][i] * target[atom][j];
            }
        }
        const double sxx = covariance[0][0], sxy = covariance[0][1], sxz = covariance[0][2];
        const double syx = covariance[1][0], syy = covariance[1][1], syz = covariance[1][2];
        const double szx = covariance[2][0], szy = covariance[2][1], szz = covariance[2][2];
        const std::array<std::array<double, 4>, 4> horn{{
            {{sxx + syy + szz, syz - szy, szx - sxz, sxy - syx}},
            {{syz - szy, sxx - syy - szz, sxy + syx, szx + sxz}},
            {{szx - sxz, sxy + syx, -sxx + syy - szz, syz + szy}},
            {{sxy - syx, szx + sxz, syz + szy, -sxx - syy + szz}},
        }};
        const auto quaternion = largest_eigenvector(horn);
        for (Vec3& point : current) point = rotate(point, quaternion);
    }
    double squared = 0.0;
    for (std::size_t atom = 0; atom < current.size(); ++atom) {
        const Vec3 delta = subtract(current[atom], target[atom]);
        squared += dot(delta, delta);
    }
    return std::sqrt(squared / static_cast<double>(current.size()));
}

double coordination_number(const Frame& frame, const std::vector<std::size_t>& centers,
                           const std::vector<std::size_t>& neighbors, double cutoff_nm) {
    validate_selection(frame, centers, "Coordination center");
    validate_selection(frame, neighbors, "Coordination neighbor");
    const auto found = contacts(frame, centers, neighbors, cutoff_nm);
    std::map<std::size_t, std::vector<std::size_t>> unique;
    for (const Contact& contact : found) {
        auto& values = unique[contact.atom_a];
        if (std::find(values.begin(), values.end(), contact.atom_b) == values.end()) values.push_back(contact.atom_b);
    }
    std::size_t total = 0;
    for (const std::size_t center : centers) total += unique[center].size();
    return static_cast<double>(total) / static_cast<double>(centers.size());
}

bool aromatic_stacking(const Frame& frame, const std::vector<std::size_t>& first_ring,
                       const std::vector<std::size_t>& second_ring, double centroid_cutoff_nm,
                       double maximum_normal_angle_degrees) {
    if (first_ring.size() < 3 || second_ring.size() < 3) {
        throw std::invalid_argument("AROMATIC_STACKING requires at least three atoms in each ring selection");
    }
    if (!std::isfinite(centroid_cutoff_nm) || centroid_cutoff_nm < 0.0 ||
        !std::isfinite(maximum_normal_angle_degrees) || maximum_normal_angle_degrees < 0.0 ||
        maximum_normal_angle_degrees > 90.0) {
        throw std::invalid_argument("Invalid aromatic stacking distance or angle criterion");
    }
    const auto first = unwrapped(frame, first_ring), second = unwrapped(frame, second_ring);
    const Vec3 first_center = centroid(first), second_center = centroid(second);
    const double separation = norm(displacement(frame, first_center, second_center));
    return separation <= centroid_cutoff_nm + 1e-12 &&
           vector_angle(ring_normal(first), ring_normal(second), true) <=
               maximum_normal_angle_degrees + 1e-12;
}

double molecular_orientation(const Frame& frame, const std::vector<std::size_t>& origin,
                             const std::vector<std::size_t>& target, char axis) {
    const Vec3 first = centroid(unwrapped(frame, origin));
    const Vec3 second = centroid(unwrapped(frame, target));
    Vec3 reference{};
    if (axis == 'x' || axis == 'X') reference[0] = 1.0;
    else if (axis == 'y' || axis == 'Y') reference[1] = 1.0;
    else if (axis == 'z' || axis == 'Z') reference[2] = 1.0;
    else throw std::invalid_argument("ORIENTATION axis must be x, y, or z");
    return vector_angle(displacement(frame, first, second), reference, true);
}

} // namespace ensembleql
