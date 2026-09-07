#include "ensembleql/observable.hpp"
#include "ensembleql/geometry.hpp"
#include "ensembleql/molecular.hpp"

#include <set>
#include <utility>

namespace ensembleql {

DistanceObservable::DistanceObservable(std::vector<std::size_t> a, std::vector<std::size_t> b)
    : a_(std::move(a)), b_(std::move(b)) {}
double DistanceObservable::evaluate(const Frame& frame) const { return minimum_distance(frame, a_, b_); }

ContactObservable::ContactObservable(std::vector<std::size_t> a, std::vector<std::size_t> b, double cutoff_nm)
    : neighbor_list_(std::move(a), std::move(b), cutoff_nm) {}
double ContactObservable::evaluate(const Frame& frame) const {
    return neighbor_list_.has_contact(frame) ? 1.0 : 0.0;
}

ContactCountObservable::ContactCountObservable(std::vector<std::size_t> a, std::vector<std::size_t> b,
                                               double cutoff_nm, ContactMode mode, Topology topology)
    : mode_(mode), topology_(std::move(topology)),
      neighbor_list_(std::move(a), std::move(b), cutoff_nm) {}
double ContactCountObservable::evaluate(const Frame& frame) const {
    const auto atom_contacts = neighbor_list_.contacts(frame);
    std::set<Bond> unique_atom_contacts;
    for (const Contact& contact : atom_contacts) {
        Bond pair{contact.atom_a, contact.atom_b};
        if (pair[1] < pair[0]) std::swap(pair[0], pair[1]);
        unique_atom_contacts.insert(pair);
    }
    if (mode_ == ContactMode::Atom) return static_cast<double>(unique_atom_contacts.size());
    using Residue = std::pair<char, int>;
    std::set<std::pair<Residue, Residue>> residue_contacts;
    for (const Bond& contact : unique_atom_contacts) {
        Residue first{topology_.atoms().at(contact[0]).chain,
                      topology_.atoms().at(contact[0]).residue_index};
        Residue second{topology_.atoms().at(contact[1]).chain,
                       topology_.atoms().at(contact[1]).residue_index};
        if (first == second) continue;
        if (second < first) std::swap(first, second);
        residue_contacts.emplace(first, second);
    }
    return static_cast<double>(residue_contacts.size());
}

RgObservable::RgObservable(std::vector<std::size_t> selection, std::vector<Bond> bonds,
                           std::vector<double> masses_da)
    : selection_(std::move(selection)), bonds_(std::move(bonds)), masses_da_(std::move(masses_da)) {}
double RgObservable::evaluate(const Frame& frame) const {
    if (!masses_da_.empty()) {
        return mass_weighted_radius_of_gyration(frame, selection_, bonds_, masses_da_);
    }
    return radius_of_gyration(frame, selection_, bonds_);
}

HydrogenBondObservable::HydrogenBondObservable(
    std::vector<std::size_t> donors, std::vector<std::size_t> acceptors,
    Topology topology, double distance_cutoff_nm, double minimum_angle_degrees)
    : donors_(std::move(donors)), acceptors_(std::move(acceptors)),
      topology_(std::move(topology)), distance_cutoff_nm_(distance_cutoff_nm),
      minimum_angle_degrees_(minimum_angle_degrees) {}
double HydrogenBondObservable::evaluate(const Frame& frame) const {
    return hydrogen_bond(frame, donors_, acceptors_, topology_, distance_cutoff_nm_,
                         minimum_angle_degrees_) ? 1.0 : 0.0;
}

DihedralObservable::DihedralObservable(std::size_t first, std::size_t second,
                                       std::size_t third, std::size_t fourth)
    : first_(first), second_(second), third_(third), fourth_(fourth) {}
double DihedralObservable::evaluate(const Frame& frame) const {
    return dihedral_degrees(frame, first_, second_, third_, fourth_);
}

HelixObservable::HelixObservable(std::vector<std::size_t> selection,
                                 Topology topology, double minimum_fraction)
    : selection_(std::move(selection)), topology_(std::move(topology)),
      minimum_fraction_(minimum_fraction) {}
double HelixObservable::evaluate(const Frame& frame) const {
    return helix_fraction(frame, selection_, topology_) >= minimum_fraction_ ? 1.0 : 0.0;
}

SasaObservable::SasaObservable(std::vector<std::size_t> selection, Topology topology,
                               double probe_radius_nm, std::size_t sphere_points)
    : selection_(std::move(selection)), topology_(std::move(topology)),
      probe_radius_nm_(probe_radius_nm), sphere_points_(sphere_points) {}
double SasaObservable::evaluate(const Frame& frame) const {
    return solvent_accessible_surface_area(frame, selection_, topology_, probe_radius_nm_,
                                           sphere_points_);
}

RmsdObservable::RmsdObservable(std::vector<std::size_t> selection, bool align)
    : selection_(std::move(selection)), align_(align) {}
double RmsdObservable::evaluate(const Frame& frame) const {
    if (reference_.empty()) {
        reference_.reserve(selection_.size());
        for (const std::size_t index : selection_) reference_.push_back(frame.coordinates.at(index));
    }
    return aligned_rmsd(frame, selection_, reference_, align_);
}

CoordinationNumberObservable::CoordinationNumberObservable(
    std::vector<std::size_t> centers, std::vector<std::size_t> neighbors,
    double cutoff_nm)
    : centers_(std::move(centers)), neighbors_(std::move(neighbors)),
      cutoff_nm_(cutoff_nm) {}
double CoordinationNumberObservable::evaluate(const Frame& frame) const {
    return coordination_number(frame, centers_, neighbors_, cutoff_nm_);
}

SaltBridgeObservable::SaltBridgeObservable(std::vector<std::size_t> positive,
                                           std::vector<std::size_t> negative,
                                           double cutoff_nm)
    : neighbor_list_(std::move(positive), std::move(negative), cutoff_nm) {}
double SaltBridgeObservable::evaluate(const Frame& frame) const {
    return neighbor_list_.has_contact(frame) ? 1.0 : 0.0;
}

AromaticStackingObservable::AromaticStackingObservable(
    std::vector<std::size_t> first_ring, std::vector<std::size_t> second_ring,
    double centroid_cutoff_nm, double maximum_angle_degrees)
    : first_ring_(std::move(first_ring)), second_ring_(std::move(second_ring)),
      centroid_cutoff_nm_(centroid_cutoff_nm),
      maximum_angle_degrees_(maximum_angle_degrees) {}
double AromaticStackingObservable::evaluate(const Frame& frame) const {
    return aromatic_stacking(frame, first_ring_, second_ring_, centroid_cutoff_nm_,
                             maximum_angle_degrees_) ? 1.0 : 0.0;
}

SurfaceDistanceObservable::SurfaceDistanceObservable(
    std::vector<std::size_t> molecule, std::vector<std::size_t> surface)
    : molecule_(std::move(molecule)), surface_(std::move(surface)) {}
double SurfaceDistanceObservable::evaluate(const Frame& frame) const {
    return minimum_distance(frame, molecule_, surface_);
}

OrientationObservable::OrientationObservable(std::vector<std::size_t> origin,
                                             std::vector<std::size_t> target,
                                             char axis)
    : origin_(std::move(origin)), target_(std::move(target)), axis_(axis) {}
double OrientationObservable::evaluate(const Frame& frame) const {
    return molecular_orientation(frame, origin_, target_, axis_);
}

} // namespace ensembleql
