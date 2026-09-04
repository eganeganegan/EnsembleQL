#pragma once

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
};

class Topology {
public:
    Topology() = default;
    explicit Topology(std::vector<Atom> atoms);
    static Topology from_pdb(const std::string& path);
    const std::vector<Atom>& atoms() const noexcept { return atoms_; }
    std::size_t size() const noexcept { return atoms_.size(); }
private:
    std::vector<Atom> atoms_;
};

} // namespace ensembleql
