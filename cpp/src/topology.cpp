#include "ensembleql/topology.hpp"

#include <algorithm>
#include <cctype>
#include <fstream>
#include <map>
#include <set>
#include <stdexcept>
#include <utility>

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

Topology::Topology(std::vector<Atom> atoms, std::vector<Bond> bonds)
    : atoms_(std::move(atoms)), bonds_(std::move(bonds)) {
    std::set<Bond> normalized;
    for (auto bond : bonds_) {
        if (bond[0] >= atoms_.size() || bond[1] >= atoms_.size()) {
            throw std::out_of_range("Topology bond index exceeds atom count");
        }
        if (bond[0] == bond[1]) throw std::invalid_argument("Topology bond cannot connect an atom to itself");
        if (bond[1] < bond[0]) std::swap(bond[0], bond[1]);
        normalized.insert(bond);
    }
    bonds_.assign(normalized.begin(), normalized.end());
}

Topology Topology::from_pdb(const std::string& path) {
    std::ifstream input(path);
    if (!input) throw std::runtime_error("Cannot open PDB topology: " + path);
    std::vector<Atom> atoms;
    std::map<int, std::size_t> serial_to_index;
    std::vector<std::array<int, 2>> serial_bonds;
    std::string line;
    while (std::getline(input, line)) {
        if (line.rfind("CONECT", 0) == 0) {
            try {
                const int source = std::stoi(field(line, 6, 5));
                for (std::size_t position = 11; position < line.size(); position += 5) {
                    const std::string target = field(line, position, 5);
                    if (!target.empty()) serial_bonds.push_back({source, std::stoi(target)});
                }
            } catch (const std::exception&) {
                throw std::runtime_error("Malformed PDB CONECT record: " + line);
            }
            continue;
        }
        if (line.rfind("ATOM  ", 0) != 0 && line.rfind("HETATM", 0) != 0) continue;
        try {
            Atom atom;
            atom.index = atoms.size();
            const int serial = std::stoi(field(line, 6, 5));
            atom.name = field(line, 12, 4);
            atom.residue_name = field(line, 17, 3);
            atom.chain = line.size() > 21 ? line[21] : ' ';
            atom.residue_index = std::stoi(field(line, 22, 4));
            atom.element = field(line, 76, 2);
            if (atom.element.empty() && !atom.name.empty()) {
                const auto it = std::find_if(atom.name.begin(), atom.name.end(), [](unsigned char c) { return std::isalpha(c); });
                if (it != atom.name.end()) atom.element = std::string(1, static_cast<char>(std::toupper(static_cast<unsigned char>(*it))));
            }
            serial_to_index[serial] = atom.index;
            atoms.push_back(std::move(atom));
        } catch (const std::exception&) {
            throw std::runtime_error("Malformed PDB atom record: " + line);
        }
    }
    if (atoms.empty()) throw std::runtime_error("PDB topology contains no ATOM/HETATM records: " + path);
    std::vector<Bond> bonds;
    for (const auto& serial_bond : serial_bonds) {
        const auto first = serial_to_index.find(serial_bond[0]);
        const auto second = serial_to_index.find(serial_bond[1]);
        // CONECT records may refer to atoms omitted from the selected model.
        if (first != serial_to_index.end() && second != serial_to_index.end() && first->second != second->second) {
            bonds.push_back({first->second, second->second});
        }
    }
    return Topology(std::move(atoms), std::move(bonds));
}

} // namespace ensembleql
