#include "ensembleql/topology.hpp"

#include <algorithm>
#include <cctype>
#include <fstream>
#include <stdexcept>

namespace ensembleql {
namespace {
std::string trim(std::string value) {
    const auto first = value.find_first_not_of(" \t\r\n");
    if (first == std::string::npos) return {};
    const auto last = value.find_last_not_of(" \t\r\n");
    return value.substr(first, last - first + 1);
}
std::string field(const std::string& line, std::size_t pos, std::size_t count) {
    return pos < line.size() ? trim(line.substr(pos, count)) : std::string{};
}
} // namespace

Topology::Topology(std::vector<Atom> atoms) : atoms_(std::move(atoms)) {}

Topology Topology::from_pdb(const std::string& path) {
    std::ifstream input(path);
    if (!input) throw std::runtime_error("Cannot open PDB topology: " + path);
    std::vector<Atom> atoms;
    std::string line;
    while (std::getline(input, line)) {
        if (line.rfind("ATOM  ", 0) != 0 && line.rfind("HETATM", 0) != 0) continue;
        try {
            Atom atom;
            atom.index = atoms.size();
            atom.name = field(line, 12, 4);
            atom.residue_name = field(line, 17, 3);
            atom.chain = line.size() > 21 ? line[21] : ' ';
            atom.residue_index = std::stoi(field(line, 22, 4));
            atom.element = field(line, 76, 2);
            if (atom.element.empty() && !atom.name.empty()) {
                const auto it = std::find_if(atom.name.begin(), atom.name.end(), [](unsigned char c) { return std::isalpha(c); });
                if (it != atom.name.end()) atom.element = std::string(1, static_cast<char>(std::toupper(static_cast<unsigned char>(*it))));
            }
            atoms.push_back(std::move(atom));
        } catch (const std::exception&) {
            throw std::runtime_error("Malformed PDB atom record: " + line);
        }
    }
    if (atoms.empty()) throw std::runtime_error("PDB topology contains no ATOM/HETATM records: " + path);
    return Topology(std::move(atoms));
}

} // namespace ensembleql
