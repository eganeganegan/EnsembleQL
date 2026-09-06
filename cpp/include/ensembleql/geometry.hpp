#pragma once

#include "ensembleql/topology.hpp"
#include "ensembleql/trajectory.hpp"

#include <cstddef>
#include <vector>

namespace ensembleql {

struct Contact {
    std::size_t atom_a{};
    std::size_t atom_b{};
    double distance_nm{};
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
