#pragma once

#include "ensembleql/topology.hpp"
#include "ensembleql/trajectory.hpp"

#include <cstddef>
#include <string>
#include <vector>

namespace ensembleql {

bool sasa_element_supported(const std::string& element);

double angle_degrees(const Frame& frame, std::size_t first,
                     std::size_t vertex, std::size_t third);
double dihedral_degrees(const Frame& frame, std::size_t first, std::size_t second,
                        std::size_t third, std::size_t fourth);
bool hydrogen_bond(const Frame& frame,
                   const std::vector<std::size_t>& donors,
                   const std::vector<std::size_t>& acceptors,
                   const Topology& topology,
                   double distance_cutoff_nm = 0.35,
                   double minimum_angle_degrees = 150.0);
double helix_fraction(const Frame& frame,
                      const std::vector<std::size_t>& selection,
                      const Topology& topology);
double solvent_accessible_surface_area(const Frame& frame,
                                       const std::vector<std::size_t>& selection,
                                       const Topology& topology,
                                       double probe_radius_nm = 0.14,
                                       std::size_t sphere_points = 96);
double aligned_rmsd(const Frame& frame,
                    const std::vector<std::size_t>& selection,
                    const std::vector<Vec3>& reference,
                    bool align = true);
double coordination_number(const Frame& frame,
                           const std::vector<std::size_t>& centers,
                           const std::vector<std::size_t>& neighbors,
                           double cutoff_nm = 0.35);
bool aromatic_stacking(const Frame& frame,
                       const std::vector<std::size_t>& first_ring,
                       const std::vector<std::size_t>& second_ring,
                       double centroid_cutoff_nm = 0.55,
                       double maximum_normal_angle_degrees = 30.0);
double molecular_orientation(const Frame& frame,
                             const std::vector<std::size_t>& origin,
                             const std::vector<std::size_t>& target,
                             char axis);

} // namespace ensembleql
