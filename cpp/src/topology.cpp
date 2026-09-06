#include "ensembleql/topology.hpp"

#include <algorithm>
#include <atomic>
#include <cctype>
#include <fstream>
#include <map>
#include <set>
#include <stdexcept>
#include <string_view>
#include <utility>

namespace ensembleql {
namespace {
std::uint64_t next_topology_identity() {
    static std::atomic<std::uint64_t> identity{1};
    return identity.fetch_add(1, std::memory_order_relaxed);
}

std::string trim(std::string value) {
    const auto first = value.find_first_not_of(" \t\r\n");
    if (first == std::string::npos) return {};
    const auto last = value.find_last_not_of(" \t\r\n");
    return value.substr(first, last - first + 1);
}
std::string field(const std::string& line, std::size_t pos, std::size_t count) {
    return pos < line.size() ? trim(line.substr(pos, count)) : std::string{};
}

std::string normalized_element(std::string element) {
    element = trim(std::move(element));
    std::transform(element.begin(), element.end(), element.begin(), [](unsigned char character) {
        return static_cast<char>(std::toupper(character));
    });
    return element;
}
} // namespace

Topology::Topology() : cache_identity_(next_topology_identity()) {}

double standard_atomic_mass_da(const std::string& element) {
    const std::string symbol = normalized_element(element);
    // CIAAW 2024 abridged standard atomic weights. Elements without a
    // standard atomic weight intentionally map to zero (unknown).
    static constexpr std::pair<std::string_view, double> masses[] = {
        {"H", 1.0080}, {"HE", 4.0026}, {"LI", 6.94}, {"BE", 9.0122},
        {"B", 10.81}, {"C", 12.011}, {"N", 14.007}, {"O", 15.999},
        {"F", 18.998}, {"NE", 20.180}, {"NA", 22.990}, {"MG", 24.305},
        {"AL", 26.982}, {"SI", 28.085}, {"P", 30.974}, {"S", 32.06},
        {"CL", 35.45}, {"AR", 39.95}, {"K", 39.098}, {"CA", 40.078},
        {"SC", 44.956}, {"TI", 47.867}, {"V", 50.942}, {"CR", 51.996},
        {"MN", 54.938}, {"FE", 55.845}, {"CO", 58.933}, {"NI", 58.693},
        {"CU", 63.546}, {"ZN", 65.38}, {"GA", 69.723}, {"GE", 72.630},
        {"AS", 74.922}, {"SE", 78.971}, {"BR", 79.904}, {"KR", 83.798},
        {"RB", 85.468}, {"SR", 87.62}, {"Y", 88.906}, {"ZR", 91.222},
        {"NB", 92.906}, {"MO", 95.95}, {"RU", 101.07}, {"RH", 102.91},
        {"PD", 106.42}, {"AG", 107.87}, {"CD", 112.41}, {"IN", 114.82},
        {"SN", 118.71}, {"SB", 121.76}, {"TE", 127.60}, {"I", 126.90},
        {"XE", 131.29}, {"CS", 132.91}, {"BA", 137.33}, {"LA", 138.91},
        {"CE", 140.12}, {"PR", 140.91}, {"ND", 144.24}, {"SM", 150.36},
        {"EU", 151.96}, {"GD", 157.25}, {"TB", 158.93}, {"DY", 162.50},
        {"HO", 164.93}, {"ER", 167.26}, {"TM", 168.93}, {"YB", 173.05},
        {"LU", 174.97}, {"HF", 178.49}, {"TA", 180.95}, {"W", 183.84},
        {"RE", 186.21}, {"OS", 190.23}, {"IR", 192.22}, {"PT", 195.08},
        {"AU", 196.97}, {"HG", 200.59}, {"TL", 204.38}, {"PB", 207.2},
        {"BI", 208.98}, {"TH", 232.04}, {"PA", 231.04}, {"U", 238.03},
    };
    for (const auto& [candidate, mass] : masses) {
        if (candidate == symbol) return mass;
    }
    return 0.0;
}

Topology::Topology(std::vector<Atom> atoms, std::vector<Bond> bonds)
    : atoms_(std::move(atoms)), bonds_(std::move(bonds)),
      cache_identity_(next_topology_identity()) {
    for (Atom& atom : atoms_) {
        if (atom.mass_da == 0.0) atom.mass_da = standard_atomic_mass_da(atom.element);
    }
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

std::vector<double> Topology::masses_da() const {
    std::vector<double> result;
    result.reserve(atoms_.size());
    for (const Atom& atom : atoms_) result.push_back(atom.mass_da);
    return result;
}

Topology Topology::from_pdb(const std::string& path) {
    std::ifstream input(path);
    if (!input) throw std::runtime_error("Cannot open PDB topology: " + path);
    std::vector<Atom> atoms;
    std::map<int, std::size_t> serial_to_index;
    std::vector<std::array<int, 2>> serial_bonds;
    std::string line;
    bool saw_model = false;
    bool first_model_active = false;
    bool first_model_complete = false;
    while (std::getline(input, line)) {
        if (line.rfind("MODEL ", 0) == 0) {
            if (!saw_model) {
                saw_model = true;
                first_model_active = true;
            } else {
                first_model_active = false;
            }
            continue;
        }
        if (line.rfind("ENDMDL", 0) == 0) {
            if (first_model_active) first_model_complete = true;
            first_model_active = false;
            continue;
        }
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
        if ((saw_model && !first_model_active) || first_model_complete) continue;
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
