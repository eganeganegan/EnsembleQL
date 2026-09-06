#include "ensembleql/engine.hpp"
#include "ensembleql/trajectory.hpp"

#include <chemfiles.hpp>

#include <cmath>
#include <filesystem>
#include <iostream>
#include <string>

using namespace ensembleql;

namespace {
int failures = 0;

void check(bool condition, const std::string& label) {
    if (!condition) {
        std::cerr << "FAIL: " << label << '\n';
        ++failures;
    }
}

bool close(double first, double second) {
    return std::abs(first - second) < 1e-5;
}

void write_trajectory(const std::string& path) {
    chemfiles::Trajectory output(path, 'w');
    for (std::size_t index = 0; index < 2; ++index) {
        chemfiles::Frame frame(chemfiles::UnitCell({20, 20, 20}));
        const double left = index == 0 ? 1.0 : 2.0;
        const double right = index == 0 ? 19.0 : 18.0;
        frame.add_atom(chemfiles::Atom("C"), {left, 0.0, 0.0});
        frame.add_atom(chemfiles::Atom("C"), {right, 0.0, 0.0});
        frame.set("time", chemfiles::Property(static_cast<double>(index)));
        output.write(frame);
    }
}

void test_format(const std::string& extension) {
    const std::string path = std::string(ENSEMBLEQL_BINARY_DIR) + "/ensembleql_reader_test" + extension;
    std::filesystem::remove(path);
    write_trajectory(path);

    ChemfilesReader reader(path, 2);
    Frame frame;
    check(reader.next(frame), extension + " first frame read");
    check(frame.coordinates.size() == 2 && close(frame.coordinates[0][0], 0.1),
          extension + " angstrom coordinate normalization (received " +
              (frame.coordinates.empty() ? std::string("no coordinates") : std::to_string(frame.coordinates[0][0])) + ")");
    check(frame.box_nm.has_value() && close((*frame.box_nm)[0], 2.0),
          extension + " orthorhombic cell normalization");
    check(close(frame.time_ps, 0.0), extension + " time property");
    check(reader.next(frame) && close(frame.time_ps, 1.0), extension + " second frame time");
    check(!reader.next(frame), extension + " streaming end");

    auto trajectory = Trajectory::from_files(
        path, std::string(ENSEMBLEQL_SOURCE_DIR) + "/tests/data/pbc.pdb");
    const auto events = Engine().query(trajectory, "FIND CONTACT(resid 1, resid 2, cutoff=0.41nm);");
    check(events.size() == 1 && close(events[0].end_time, 1.0), extension + " query integration");
    std::filesystem::remove(path);
}
} // namespace

int main() {
    check(chemfiles_backend_available(), "chemfiles backend availability");
    test_format(".dcd");
    test_format(".xtc");
    test_format(".trr");
    if (failures != 0) {
        std::cerr << failures << " chemfiles test(s) failed\n";
        return 1;
    }
    std::cout << "All EnsembleQL chemfiles tests passed\n";
    return 0;
}
