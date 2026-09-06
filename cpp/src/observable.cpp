#include "ensembleql/observable.hpp"
#include "ensembleql/geometry.hpp"

#include <algorithm>
#include <cmath>
#include <utility>

namespace ensembleql {

DistanceObservable::DistanceObservable(std::vector<std::size_t> a, std::vector<std::size_t> b)
    : a_(std::move(a)), b_(std::move(b)) {}
double DistanceObservable::evaluate(const Frame& frame) const { return minimum_distance(frame, a_, b_); }

ContactObservable::ContactObservable(std::vector<std::size_t> a, std::vector<std::size_t> b, double cutoff_nm)
    : a_(std::move(a)), b_(std::move(b)), cutoff_nm_(cutoff_nm) {}
double ContactObservable::evaluate(const Frame& frame) const {
    const double tolerance = 1e-12 * std::max(1.0, std::abs(cutoff_nm_));
    return minimum_distance(frame, a_, b_) <= cutoff_nm_ + tolerance ? 1.0 : 0.0;
}

ContactCountObservable::ContactCountObservable(std::vector<std::size_t> a, std::vector<std::size_t> b, double cutoff_nm)
    : a_(std::move(a)), b_(std::move(b)), cutoff_nm_(cutoff_nm) {}
double ContactCountObservable::evaluate(const Frame& frame) const { return static_cast<double>(contacts(frame, a_, b_, cutoff_nm_).size()); }

RgObservable::RgObservable(std::vector<std::size_t> selection, std::vector<Bond> bonds)
    : selection_(std::move(selection)), bonds_(std::move(bonds)) {}
double RgObservable::evaluate(const Frame& frame) const { return radius_of_gyration(frame, selection_, bonds_); }

} // namespace ensembleql
