#include "ensembleql/trajectory.hpp"

#include <algorithm>
#include <cmath>
#include <numbers>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace ensembleql {
namespace {
std::string trim(std::string value) {
    const auto first = value.find_first_not_of(" \t\r\n");
    if (first == std::string::npos) return {};
    const auto last = value.find_last_not_of(" \t\r\n");
    return value.substr(first, last - first + 1);
}

std::string field(const std::string& line, std::size_t position, std::size_t count) {
    return position < line.size() ? trim(line.substr(position, count)) : std::string{};
}

bool atom_record(const std::string& line) {
    return line.rfind("ATOM  ", 0) == 0 || line.rfind("HETATM", 0) == 0;
}

struct ParsedCell {
    std::optional<Vec3> box_nm;
    std::optional<PeriodicCell> cell_nm;
};

ParsedCell parse_cryst1(const std::string& line) {
    try {
        const Vec3 lengths_nm{
            std::stod(field(line, 6, 9)) * 0.1,
            std::stod(field(line, 15, 9)) * 0.1,
            std::stod(field(line, 24, 9)) * 0.1,
        };
        const Vec3 angles_deg{
            std::stod(field(line, 33, 7)),
            std::stod(field(line, 40, 7)),
            std::stod(field(line, 47, 7)),
        };
        for (const double length : lengths_nm) {
            if (!std::isfinite(length) || length <= 0.0) {
                throw std::runtime_error("cell lengths must be finite and positive");
            }
        }
        for (const double angle : angles_deg) {
            if (!std::isfinite(angle) || angle <= 0.0 || angle >= 180.0) {
                throw std::runtime_error("cell angles must be between 0 and 180 degrees");
            }
        }

        constexpr double orthogonal_tolerance = 1e-8;
        const bool orthorhombic = std::all_of(
            angles_deg.begin(), angles_deg.end(), [](double angle) {
                return std::abs(angle - 90.0) <= orthogonal_tolerance;
            });
        const bool unitary_placeholder = orthorhombic && std::all_of(
            lengths_nm.begin(), lengths_nm.end(), [](double length) {
                return std::abs(length - 0.1) <= 1e-12;
            });
        if (unitary_placeholder) return {};
        if (orthorhombic) {
            PeriodicCell cell;
            cell.vectors_nm = {
                Vec3{lengths_nm[0], 0.0, 0.0},
                Vec3{0.0, lengths_nm[1], 0.0},
                Vec3{0.0, 0.0, lengths_nm[2]},
            };
            return {lengths_nm, cell};
        }

        const double alpha = angles_deg[0] * std::numbers::pi / 180.0;
        const double beta = angles_deg[1] * std::numbers::pi / 180.0;
        const double gamma = angles_deg[2] * std::numbers::pi / 180.0;
        const double sin_gamma = std::sin(gamma);
        if (std::abs(sin_gamma) <= 1e-12) throw std::runtime_error("cell gamma angle is singular");

        const Vec3 a{lengths_nm[0], 0.0, 0.0};
        const Vec3 b{lengths_nm[1] * std::cos(gamma), lengths_nm[1] * sin_gamma, 0.0};
        const double c_x = lengths_nm[2] * std::cos(beta);
        const double c_y = lengths_nm[2] *
            (std::cos(alpha) - std::cos(beta) * std::cos(gamma)) / sin_gamma;
        double c_z_squared = lengths_nm[2] * lengths_nm[2] - c_x * c_x - c_y * c_y;
        if (c_z_squared < -1e-12) throw std::runtime_error("cell angles do not form a valid volume");
        c_z_squared = std::max(0.0, c_z_squared);
        const Vec3 c{c_x, c_y, std::sqrt(c_z_squared)};
        if (c[2] <= 1e-12) throw std::runtime_error("cell matrix is singular");
        PeriodicCell cell;
        cell.vectors_nm = {a, b, c};
        return {std::nullopt, cell};
    } catch (const std::exception& error) {
        throw std::runtime_error("Malformed PDB CRYST1 record: " + line + " (" + error.what() + ")");
    }
}

struct ParsedAtom {
    std::string identity;
    Vec3 coordinate;
};

ParsedAtom parse_atom(const std::string& line) {
    try {
        Vec3 coordinate{
            std::stod(field(line, 30, 8)) * 0.1,
            std::stod(field(line, 38, 8)) * 0.1,
            std::stod(field(line, 46, 8)) * 0.1,
        };
        for (const double value : coordinate) {
            if (!std::isfinite(value)) throw std::runtime_error("coordinates must be finite");
        }
        const std::string identity = field(line, 0, 6) + "|" + field(line, 12, 4) + "|" +
            field(line, 16, 1) + "|" + field(line, 17, 3) + "|" + field(line, 21, 1) + "|" +
            field(line, 22, 4) + "|" + field(line, 26, 1) + "|" + field(line, 76, 2);
        return {identity, coordinate};
    } catch (const std::exception& error) {
        throw std::runtime_error("Malformed PDB coordinate record: " + line + " (" + error.what() + ")");
    }
}
} // namespace

PDBReader::PDBReader(std::string path, std::size_t expected_atoms, double default_step_ps)
    : path_(std::move(path)), expected_atoms_(expected_atoms), default_step_ps_(default_step_ps) {
    if (!std::isfinite(default_step_ps_) || default_step_ps_ <= 0.0) {
        throw std::invalid_argument("PDB default timestep must be finite and positive");
    }
    reset();
}

void PDBReader::reset() {
    stream_.close();
    stream_.clear();
    stream_.open(path_);
    if (!stream_) throw std::runtime_error("Cannot open PDB trajectory: " + path_);
    frame_index_ = 0;
    saw_model_records_ = false;
    finished_single_frame_ = false;
    reference_atom_identities_.clear();
    box_nm_.reset();
    cell_nm_.reset();
}

bool PDBReader::finish_frame(Frame& frame, std::vector<Vec3> coordinates,
                             const std::vector<std::string>& atom_identities) {
    if (coordinates.size() != expected_atoms_) {
        throw std::runtime_error("PDB frame has " + std::to_string(coordinates.size()) +
                                 " atoms, topology has " + std::to_string(expected_atoms_));
    }
    if (reference_atom_identities_.empty()) {
        reference_atom_identities_ = atom_identities;
    } else if (atom_identities != reference_atom_identities_) {
        throw std::runtime_error("PDB MODEL atom identities or ordering changed between frames");
    }
    frame.time_ps = static_cast<double>(frame_index_) * default_step_ps_;
    frame.coordinates = std::move(coordinates);
    frame.box_nm = box_nm_;
    frame.cell_nm = cell_nm_;
    ++frame_index_;
    return true;
}

bool PDBReader::next(Frame& frame) {
    if (finished_single_frame_) return false;
    std::vector<Vec3> coordinates;
    std::vector<std::string> atom_identities;
    bool inside_model = false;
    bool model_started = false;
    std::string line;
    while (std::getline(stream_, line)) {
        if (line.rfind("CRYST1", 0) == 0) {
            const ParsedCell parsed = parse_cryst1(line);
            box_nm_ = parsed.box_nm;
            cell_nm_ = parsed.cell_nm;
            continue;
        }
        if (line.rfind("MODEL ", 0) == 0) {
            if (!coordinates.empty() || model_started) {
                throw std::runtime_error("Malformed PDB trajectory: nested MODEL or atoms before MODEL");
            }
            saw_model_records_ = true;
            inside_model = true;
            model_started = true;
            continue;
        }
        if (line.rfind("ENDMDL", 0) == 0) {
            if (!inside_model || coordinates.empty()) {
                throw std::runtime_error("Malformed PDB trajectory: empty or unmatched ENDMDL");
            }
            return finish_frame(frame, std::move(coordinates), atom_identities);
        }
        if (atom_record(line)) {
            if (saw_model_records_ && !inside_model) continue;
            const ParsedAtom atom = parse_atom(line);
            atom_identities.push_back(atom.identity);
            coordinates.push_back(atom.coordinate);
            continue;
        }
        if (line.rfind("END", 0) == 0) {
            if (inside_model) {
                throw std::runtime_error("Malformed PDB trajectory: MODEL is missing ENDMDL");
            }
            finished_single_frame_ = true;
            if (!coordinates.empty()) return finish_frame(frame, std::move(coordinates), atom_identities);
            return false;
        }
    }
    if (inside_model) throw std::runtime_error("Malformed PDB trajectory: MODEL is missing ENDMDL");
    finished_single_frame_ = !saw_model_records_;
    if (!coordinates.empty()) return finish_frame(frame, std::move(coordinates), atom_identities);
    return false;
}

} // namespace ensembleql
