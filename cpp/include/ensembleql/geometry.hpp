#pragma once

#include "ensembleql/topology.hpp"
#include "ensembleql/trajectory.hpp"

#include <cstddef>
#include <optional>
#include <utility>
#include <vector>

namespace ensembleql {

struct Contact {
    std::size_t atom_a{};
    std::size_t atom_b{};
    double distance_nm{};
};

// A Verlet-style candidate list for repeated contact queries over successive
// frames. Candidate pairs are built at cutoff + skin and reused until any
// selected atom moves more than half the skin or the periodic cell changes.
class NeighborList {
public:
    NeighborList(std::vector<std::size_t> selection_a,
                 std::vector<std::size_t> selection_b,
                 double cutoff_nm, double skin_nm = 0.1);

    std::vector<Contact> contacts(const Frame& frame);
    bool has_contact(const Frame& frame);
    std::size_t rebuild_count() const noexcept { return rebuild_count_; }
    double cutoff_nm() const noexcept { return cutoff_nm_; }
    double skin_nm() const noexcept { return skin_nm_; }

private:
    void update(const Frame& frame);
    bool needs_rebuild(const Frame& frame) const;
    void rebuild(const Frame& frame);

    std::vector<std::size_t> selection_a_, selection_b_, tracked_atoms_;
    std::vector<std::pair<std::size_t, std::size_t>> candidate_pairs_;
    std::vector<Vec3> reference_coordinates_;
    std::optional<PeriodicCell> reference_cell_;
    double cutoff_nm_{};
    double skin_nm_{};
    std::size_t reference_coordinate_count_{};
    std::size_t rebuild_count_{};
    bool initialized_{};
};

double distance(const Vec3& a, const Vec3& b);
Vec3 minimum_image_displacement(const Vec3& displacement, const PeriodicCell& cell);
double minimum_image_distance(const Vec3& a, const Vec3& b, const Vec3& box_nm);
double minimum_image_distance_cell(const Vec3& a, const Vec3& b, const PeriodicCell& cell);
double minimum_distance(const Frame& frame,
                        const std::vector<std::size_t>& selection_a,
                        const std::vector<std::size_t>& selection_b);
std::vector<Contact> contacts(const Frame& frame,
                              const std::vector<std::size_t>& selection_a,
                              const std::vector<std::size_t>& selection_b,
                              double cutoff_nm);
bool has_contact(const Frame& frame,
                 const std::vector<std::size_t>& selection_a,
                 const std::vector<std::size_t>& selection_b,
                 double cutoff_nm);
std::size_t contact_count(const Frame& frame,
                          const std::vector<std::size_t>& selection_a,
                          const std::vector<std::size_t>& selection_b,
                          double cutoff_nm,
                          ContactMode mode,
                          const Topology& topology);
double radius_of_gyration(const Frame& frame,
                          const std::vector<std::size_t>& selection);
double radius_of_gyration(const Frame& frame,
                          const std::vector<std::size_t>& selection,
                          const std::vector<Bond>& bonds);
double mass_weighted_radius_of_gyration(const Frame& frame,
                                        const std::vector<std::size_t>& selection,
                                        const std::vector<Bond>& bonds,
                                        const std::vector<double>& masses_da);

} // namespace ensembleql
