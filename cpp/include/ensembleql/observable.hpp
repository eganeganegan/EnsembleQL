#pragma once

#include "ensembleql/topology.hpp"
#include "ensembleql/trajectory.hpp"

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
    std::vector<std::size_t> a_, b_;
    double cutoff_nm_;
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
    std::vector<std::size_t> a_, b_;
    double cutoff_nm_;
    ContactMode mode_;
    Topology topology_;
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

} // namespace ensembleql
