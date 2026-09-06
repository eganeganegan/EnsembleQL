#include "ensembleql/trajectory.hpp"

#include <chemfiles.hpp>

#include <memory>
#include <stdexcept>
#include <string>
#include <utility>

namespace ensembleql {

struct ChemfilesReader::Impl {
    Impl(std::string input_path, std::size_t atoms, double step_ps)
        : path(std::move(input_path)), expected_atoms(atoms), default_step_ps(step_ps) {}

    std::string path;
    std::size_t expected_atoms;
    double default_step_ps;
    std::size_t frame_index{};
    std::unique_ptr<chemfiles::Trajectory> trajectory;
};

ChemfilesReader::ChemfilesReader(std::string path, std::size_t expected_atoms, double default_step_ps)
    : impl_(std::make_unique<Impl>(std::move(path), expected_atoms, default_step_ps)) {
    reset();
}

ChemfilesReader::~ChemfilesReader() = default;
ChemfilesReader::ChemfilesReader(ChemfilesReader&&) noexcept = default;
ChemfilesReader& ChemfilesReader::operator=(ChemfilesReader&&) noexcept = default;

void ChemfilesReader::reset() {
    impl_->trajectory = std::make_unique<chemfiles::Trajectory>(impl_->path);
    impl_->frame_index = 0;
}

bool ChemfilesReader::next(Frame& frame) {
    if (impl_->trajectory->done()) return false;
    const chemfiles::Frame source = impl_->trajectory->read();
    if (source.size() != impl_->expected_atoms) {
        throw std::runtime_error("Trajectory frame has " + std::to_string(source.size()) +
                                 " atoms, topology has " + std::to_string(impl_->expected_atoms));
    }

    frame.time_ps = static_cast<double>(impl_->frame_index) * impl_->default_step_ps;
    if (const auto time = source.get("time")) frame.time_ps = time->as_double();

    frame.coordinates.clear();
    frame.coordinates.reserve(source.size());
    for (const auto& position : source.positions()) {
        frame.coordinates.push_back({position[0] * 0.1, position[1] * 0.1, position[2] * 0.1});
    }

    const auto& cell = source.cell();
    if (cell.shape() == chemfiles::UnitCell::INFINITE) {
        frame.box_nm.reset();
    } else if (cell.shape() == chemfiles::UnitCell::ORTHORHOMBIC) {
        const auto lengths = cell.lengths();
        frame.box_nm = Vec3{lengths[0] * 0.1, lengths[1] * 0.1, lengths[2] * 0.1};
    } else {
        throw std::runtime_error("Triclinic trajectory boxes are not yet supported by EnsembleQL");
    }

    ++impl_->frame_index;
    return true;
}

} // namespace ensembleql
