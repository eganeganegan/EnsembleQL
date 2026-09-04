#pragma once

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
double minimum_distance(const Frame& frame,
                        const std::vector<std::size_t>& selection_a,
                        const std::vector<std::size_t>& selection_b);
std::vector<Contact> contacts(const Frame& frame,
                              const std::vector<std::size_t>& selection_a,
                              const std::vector<std::size_t>& selection_b,
                              double cutoff_nm);
double radius_of_gyration(const Frame& frame,
                          const std::vector<std::size_t>& selection);

} // namespace ensembleql
