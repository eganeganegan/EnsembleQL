#pragma once

#include "ensembleql/topology.hpp"
#include "ensembleql/trajectory.hpp"
#include "ensembleql/geometry.hpp"

#include <cstddef>
#include <memory>
#include <string>
#include <vector>

namespace ensembleql {

class Observable {
public:
    virtual ~Observable() = default;
    virtual double evaluate(const Frame& frame) const = 0;
    virtual std::string name() const = 0;
};

class DistanceObservable final : public Observable {
public:
    DistanceObservable(std::vector<std::size_t> a, std::vector<std::size_t> b);
    double evaluate(const Frame& frame) const override;
    std::string name() const override { return "DISTANCE"; }
private:
    std::vector<std::size_t> a_, b_;
};

class ContactObservable final : public Observable {
public:
    ContactObservable(std::vector<std::size_t> a, std::vector<std::size_t> b,
                      double cutoff_nm = 0.45);
    double evaluate(const Frame& frame) const override;
    std::string name() const override { return "CONTACT"; }
private:
    mutable NeighborList neighbor_list_;
};

class ContactCountObservable final : public Observable {
public:
    ContactCountObservable(std::vector<std::size_t> a, std::vector<std::size_t> b,
                           double cutoff_nm = 0.45,
                           ContactMode mode = ContactMode::Atom,
                           Topology topology = {});
    double evaluate(const Frame& frame) const override;
    std::string name() const override { return "CONTACT_COUNT"; }
private:
    ContactMode mode_;
    Topology topology_;
    mutable NeighborList neighbor_list_;
};

class RgObservable final : public Observable {
public:
    explicit RgObservable(std::vector<std::size_t> selection, std::vector<Bond> bonds = {},
                          std::vector<double> masses_da = {});
    double evaluate(const Frame& frame) const override;
    std::string name() const override { return "RG"; }
private:
    std::vector<std::size_t> selection_;
    std::vector<Bond> bonds_;
    std::vector<double> masses_da_;
};

class HydrogenBondObservable final : public Observable {
public:
    HydrogenBondObservable(std::vector<std::size_t> donors,
                           std::vector<std::size_t> acceptors, Topology topology,
                           double distance_cutoff_nm = 0.35,
                           double minimum_angle_degrees = 150.0);
    double evaluate(const Frame& frame) const override;
    std::string name() const override { return "HBOND"; }
private:
    std::vector<std::size_t> donors_, acceptors_;
    Topology topology_;
    double distance_cutoff_nm_, minimum_angle_degrees_;
};

class DihedralObservable final : public Observable {
public:
    DihedralObservable(std::size_t first, std::size_t second,
                       std::size_t third, std::size_t fourth);
    double evaluate(const Frame& frame) const override;
    std::string name() const override { return "DIHEDRAL"; }
private:
    std::size_t first_, second_, third_, fourth_;
};

class HelixObservable final : public Observable {
public:
    HelixObservable(std::vector<std::size_t> selection, Topology topology,
                    double minimum_fraction = 0.5);
    double evaluate(const Frame& frame) const override;
    std::string name() const override { return "HELIX"; }
private:
    std::vector<std::size_t> selection_;
    Topology topology_;
    double minimum_fraction_;
};

class SasaObservable final : public Observable {
public:
    SasaObservable(std::vector<std::size_t> selection, Topology topology,
                   double probe_radius_nm = 0.14, std::size_t sphere_points = 96);
    double evaluate(const Frame& frame) const override;
    std::string name() const override { return "SASA"; }
private:
    std::vector<std::size_t> selection_;
    Topology topology_;
    double probe_radius_nm_;
    std::size_t sphere_points_;
};

class RmsdObservable final : public Observable {
public:
    RmsdObservable(std::vector<std::size_t> selection, bool align = true);
    double evaluate(const Frame& frame) const override;
    std::string name() const override { return "RMSD"; }
private:
    std::vector<std::size_t> selection_;
    bool align_;
    mutable std::vector<Vec3> reference_;
};

class CoordinationNumberObservable final : public Observable {
public:
    CoordinationNumberObservable(std::vector<std::size_t> centers,
                                 std::vector<std::size_t> neighbors,
                                 double cutoff_nm = 0.35);
    double evaluate(const Frame& frame) const override;
    std::string name() const override { return "COORDINATION_NUMBER"; }
private:
    std::vector<std::size_t> centers_, neighbors_;
    double cutoff_nm_;
};

class SaltBridgeObservable final : public Observable {
public:
    SaltBridgeObservable(std::vector<std::size_t> positive,
                         std::vector<std::size_t> negative,
                         double cutoff_nm = 0.4);
    double evaluate(const Frame& frame) const override;
    std::string name() const override { return "SALT_BRIDGE"; }
private:
    mutable NeighborList neighbor_list_;
};

class AromaticStackingObservable final : public Observable {
public:
    AromaticStackingObservable(std::vector<std::size_t> first_ring,
                              std::vector<std::size_t> second_ring,
                              double centroid_cutoff_nm = 0.55,
                              double maximum_angle_degrees = 30.0);
    double evaluate(const Frame& frame) const override;
    std::string name() const override { return "AROMATIC_STACKING"; }
private:
    std::vector<std::size_t> first_ring_, second_ring_;
    double centroid_cutoff_nm_, maximum_angle_degrees_;
};

class SurfaceDistanceObservable final : public Observable {
public:
    SurfaceDistanceObservable(std::vector<std::size_t> molecule,
                              std::vector<std::size_t> surface);
    double evaluate(const Frame& frame) const override;
    std::string name() const override { return "SURFACE_DISTANCE"; }
private:
    std::vector<std::size_t> molecule_, surface_;
};

class OrientationObservable final : public Observable {
public:
    OrientationObservable(std::vector<std::size_t> origin,
                          std::vector<std::size_t> target, char axis);
    double evaluate(const Frame& frame) const override;
    std::string name() const override { return "ORIENTATION"; }
private:
    std::vector<std::size_t> origin_, target_;
    char axis_;
};

} // namespace ensembleql
