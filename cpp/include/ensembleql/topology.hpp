#pragma once

#include <array>
#include <cstddef>
#include <string>
#include <vector>

namespace ensembleql {

struct Atom {
    std::size_t index{};
    std::string name;
    std::string residue_name;
    int residue_index{};
    char chain{' '};
    std::string element;
    double mass_da{};
};
using Bond = std::array<std::size_t, 2>;

enum class ContactMode { Atom, Residue };

class Topology {
public:
    Topology() = default;
    explicit Topology(std::vector<Atom> atoms, std::vector<Bond> bonds = {});
    static Topology from_pdb(const std::string& path);
    const std::vector<Atom>& atoms() const noexcept { return atoms_; }
    const std::vector<Bond>& bonds() const noexcept { return bonds_; }
    std::vector<double> masses_da() const;
    std::size_t size() const noexcept { return atoms_.size(); }
private:
    std::vector<Atom> atoms_;
    std::vector<Bond> bonds_;
};

double standard_atomic_mass_da(const std::string& element);

} // namespace ensembleql
