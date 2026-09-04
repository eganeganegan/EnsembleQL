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

struct Frame {
    double time_ps{};
    std::vector<Vec3> coordinates;
    std::optional<Vec3> box_nm;
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

class Trajectory {
public:
    Trajectory(Topology topology, std::shared_ptr<FrameReader> reader);
    static Trajectory from_files(const std::string& trajectory_path,
                                 const std::string& topology_path);
    const Topology& topology() const noexcept { return topology_; }
    FrameReader& reader() noexcept { return *reader_; }
private:
    Topology topology_;
    std::shared_ptr<FrameReader> reader_;
};

} // namespace ensembleql
