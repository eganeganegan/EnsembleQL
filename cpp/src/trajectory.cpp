#include "ensembleql/trajectory.hpp"
#include "ensembleql/units.hpp"

#include <filesystem>
#include <cmath>
#include <optional>
#include <regex>
#include <sstream>
#include <stdexcept>

namespace ensembleql {
namespace {
void validate_box(const Vec3& box) {
    for (const double length : box) {
        if (!std::isfinite(length) || length <= 0.0) {
            throw std::runtime_error("XYZ orthorhombic box lengths must be finite and positive");
        }
    }
}

std::optional<Vec3> parse_box(const std::string& comment) {
    const std::regex lattice_re(R"eql(Lattice\s*=\s*"([^"]+)")eql", std::regex::icase);
    std::smatch match;
    if (std::regex_search(comment, match, lattice_re)) {
        std::istringstream values(match[1].str());
        std::array<double, 9> matrix{};
        for (double& value : matrix) {
            if (!(values >> value)) throw std::runtime_error("Malformed extended XYZ Lattice matrix");
        }
        double extra{};
        if (values >> extra) throw std::runtime_error("Extended XYZ Lattice must contain exactly nine values");
        constexpr double tolerance = 1e-12;
        for (const std::size_t index : {1U, 2U, 3U, 5U, 6U, 7U}) {
            if (std::abs(matrix[index]) > tolerance) {
                throw std::runtime_error("Triclinic XYZ boxes are not yet supported; expected an orthorhombic Lattice matrix");
            }
        }
        Vec3 box{matrix[0] * 0.1, matrix[4] * 0.1, matrix[8] * 0.1};
        validate_box(box);
        return box;
    }
    const std::regex lattice_marker(R"(\bLattice\s*=)", std::regex::icase);
    if (std::regex_search(comment, lattice_marker)) {
        throw std::runtime_error("Malformed extended XYZ Lattice metadata");
    }

    const std::regex box_re(
        R"eql(\bbox\s*=\s*([+\-0-9.eE]+)\s*[,x]\s*([+\-0-9.eE]+)\s*[,x]\s*([+\-0-9.eE]+)\s*(nm|angstrom|A)\b)eql",
        std::regex::icase);
    if (std::regex_search(comment, match, box_re)) {
        Vec3 box{};
        for (std::size_t index = 0; index < 3; ++index) {
            box[index] = parse_distance(match[index + 1].str() + match[4].str()).nm;
        }
        validate_box(box);
        return box;
    }
    const std::regex box_marker(R"(\bbox\s*=)", std::regex::icase);
    if (std::regex_search(comment, box_marker)) {
        throw std::runtime_error("Malformed XYZ box metadata; expected box=Lx,Ly,Lz<unit>");
    }
    return std::nullopt;
}
} // namespace

XYZReader::XYZReader(std::string path, std::size_t expected_atoms, double default_step_ps)
    : path_(std::move(path)), expected_atoms_(expected_atoms), default_step_ps_(default_step_ps) { reset(); }

void XYZReader::reset() {
    stream_.close();
    stream_.clear();
    stream_.open(path_);
    if (!stream_) throw std::runtime_error("Cannot open XYZ trajectory: " + path_);
    frame_index_ = 0;
}

bool XYZReader::next(Frame& frame) {
    std::string line;
    while (std::getline(stream_, line) && line.empty()) {}
    if (!stream_) return false;
    std::size_t atom_count{};
    try { atom_count = static_cast<std::size_t>(std::stoul(line)); }
    catch (const std::exception&) { throw std::runtime_error("Malformed XYZ atom count: " + line); }
    if (atom_count != expected_atoms_) {
        throw std::runtime_error("XYZ frame has " + std::to_string(atom_count) +
                                 " atoms, topology has " + std::to_string(expected_atoms_));
    }
    std::string comment;
    if (!std::getline(stream_, comment)) throw std::runtime_error("XYZ frame is missing its comment line");
    frame.time_ps = static_cast<double>(frame_index_) * default_step_ps_;
    const std::regex time_re(R"((?:time\s*=\s*)?([0-9]+(?:\.[0-9]+)?)\s*(fs|ps|ns|us)\b)", std::regex::icase);
    std::smatch match;
    if (std::regex_search(comment, match, time_re)) frame.time_ps = parse_duration(match[1].str() + match[2].str()).ps;
    frame.box_nm = parse_box(comment);
    frame.coordinates.clear();
    frame.coordinates.reserve(atom_count);
    for (std::size_t i = 0; i < atom_count; ++i) {
        if (!std::getline(stream_, line)) throw std::runtime_error("Truncated XYZ frame");
        std::istringstream row(line);
        std::string symbol;
        Vec3 coordinate{};
        if (!(row >> symbol >> coordinate[0] >> coordinate[1] >> coordinate[2])) {
            throw std::runtime_error("Malformed XYZ coordinate: " + line);
        }
        for (double& value : coordinate) value *= 0.1; // XYZ convention: angstrom -> nm
        frame.coordinates.push_back(coordinate);
    }
    ++frame_index_;
    return true;
}

Trajectory::Trajectory(Topology topology, std::shared_ptr<FrameReader> reader)
    : topology_(std::move(topology)), reader_(std::move(reader)) {
    if (!reader_) throw std::invalid_argument("Trajectory requires a frame reader");
}

Trajectory Trajectory::from_files(const std::string& trajectory_path, const std::string& topology_path) {
    Topology topology = Topology::from_pdb(topology_path);
    std::string extension = std::filesystem::path(trajectory_path).extension().string();
    if (extension != ".xyz" && extension != ".XYZ") {
        throw std::runtime_error("Unsupported trajectory format '" + extension + "'; currently supported: XYZ");
    }
    auto reader = std::make_shared<XYZReader>(trajectory_path, topology.size());
    return Trajectory(std::move(topology), std::move(reader));
}

} // namespace ensembleql
