#pragma once

#include "ensembleql/topology.hpp"

#include <array>
#include <cstddef>
#include <fstream>
#include <memory>
#include <optional>
#include <string>
#include <vector>

namespace ensembleql {

using Vec3 = std::array<double, 3>;

struct PeriodicCell {
    // Cell vectors a, b, and c in nm. Cartesian positions are H * fractional,
    // where these vectors are the columns of H.
    std::array<Vec3, 3> vectors_nm{};
};

struct Frame {
    double time_ps{};
    std::vector<Vec3> coordinates;
    // Retained for source compatibility with callers constructing
    // orthorhombic frames directly. New readers also populate cell_nm.
    std::optional<Vec3> box_nm;
    std::optional<PeriodicCell> cell_nm;
};

class FrameReader {
public:
    virtual ~FrameReader() = default;
    virtual bool next(Frame& frame) = 0;
    virtual void reset() = 0;
};

class XYZReader final : public FrameReader {
public:
    XYZReader(std::string path, std::size_t expected_atoms, double default_step_ps = 1.0);
    bool next(Frame& frame) override;
    void reset() override;
private:
    std::string path_;
    std::size_t expected_atoms_;
    double default_step_ps_;
    std::size_t frame_index_{};
    std::ifstream stream_;
};

class PDBReader final : public FrameReader {
public:
    PDBReader(std::string path, std::size_t expected_atoms, double default_step_ps = 1.0);
    bool next(Frame& frame) override;
    void reset() override;
private:
    bool finish_frame(Frame& frame, std::vector<Vec3> coordinates,
                      const std::vector<std::string>& atom_identities);

    std::string path_;
    std::size_t expected_atoms_;
    double default_step_ps_;
    std::size_t frame_index_{};
    std::ifstream stream_;
    bool saw_model_records_{};
    bool finished_single_frame_{};
    std::vector<std::string> reference_atom_identities_;
    std::optional<Vec3> box_nm_;
    std::optional<PeriodicCell> cell_nm_;
};

#ifdef ENSEMBLEQL_HAS_CHEMFILES
class ChemfilesReader final : public FrameReader {
public:
    ChemfilesReader(std::string path, std::size_t expected_atoms, double default_step_ps = 1.0);
    ~ChemfilesReader() override;
    ChemfilesReader(const ChemfilesReader&) = delete;
    ChemfilesReader& operator=(const ChemfilesReader&) = delete;
    ChemfilesReader(ChemfilesReader&&) noexcept;
    ChemfilesReader& operator=(ChemfilesReader&&) noexcept;
    bool next(Frame& frame) override;
    void reset() override;
private:
    struct Impl;
    std::unique_ptr<Impl> impl_;
};
#endif

bool chemfiles_backend_available() noexcept;

class Trajectory {
public:
    Trajectory(Topology topology, std::shared_ptr<FrameReader> reader);
    static Trajectory from_files(const std::string& trajectory_path,
                                 const std::string& topology_path,
                                 double default_timestep_ps = 1.0);
    const Topology& topology() const noexcept { return topology_; }
    FrameReader& reader() noexcept { return *reader_; }
private:
    Topology topology_;
    std::shared_ptr<FrameReader> reader_;
};

} // namespace ensembleql
