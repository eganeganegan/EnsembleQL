#include "ensembleql/engine.hpp"
#include "ensembleql/event.hpp"
#include "ensembleql/geometry.hpp"
#include "ensembleql/parser.hpp"
#include "ensembleql/planner.hpp"
#include "ensembleql/selection.hpp"
#include "ensembleql/temporal.hpp"
#include "ensembleql/topology.hpp"
#include "ensembleql/trajectory.hpp"
#include "ensembleql/units.hpp"

#include <cmath>
#include <functional>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

using namespace ensembleql;

namespace {
int failures = 0;
void check(bool condition, const std::string& label) {
    if (!condition) { std::cerr << "FAIL: " << label << '\n'; ++failures; }
}
template <typename Exception = std::exception>
void check_throws(const std::function<void()>& action, const std::string& label) {
    try { action(); std::cerr << "FAIL: " << label << " (did not throw)\n"; ++failures; }
    catch (const Exception&) {}
}
bool close(double a, double b) { return std::abs(a - b) < 1e-9; }

Topology simple_topology() {
    return Topology({
        {0, "CA", "ARG", 17, 'A', "C"},
        {1, "CB", "ARG", 17, 'A', "C"},
        {2, "CA", "ASP", 42, 'A', "C"},
        {3, "CA", "GLU", 53, 'B', "C"},
        {4, "O", "HOH", 99, 'W', "O"},
    });
}

class VectorReader final : public FrameReader {
public:
    explicit VectorReader(std::vector<Frame> frames) : frames_(std::move(frames)) {}
    bool next(Frame& frame) override {
        if (position_ == frames_.size()) return false;
        frame = frames_[position_++];
        return true;
    }
    void reset() override { position_ = 0; }
private:
    std::vector<Frame> frames_;
    std::size_t position_{};
};

void test_units() {
    check(close(parse_distance("10A").nm, 1.0), "angstrom conversion");
    check(close(parse_distance("2 nm").nm, 2.0), "nm conversion");
    check(close(parse_duration("5ns").ps, 5000.0), "ns conversion");
    check(close(parse_duration("10fs").ps, 0.01), "fs conversion");
    check(close(parse_duration("2NS").ps, 2000.0), "case-insensitive time unit");
    check(close(parse_distance("2a").nm, 0.2), "case-insensitive angstrom unit");
    check_throws([&] { (void)parse_distance("2ns"); }, "distance rejects time unit");
    check_throws([&] { (void)parse_duration("2nm"); }, "time rejects distance unit");
}

void test_io_and_selections() {
    const std::string root = ENSEMBLEQL_SOURCE_DIR;
    auto topology = Topology::from_pdb(root + "/examples/idr_contact_switching/switching.pdb");
    check(topology.size() == 3, "PDB atom count");
    check(topology.atoms()[0].residue_index == 17 && topology.atoms()[1].residue_name == "ASP", "PDB fields");
    XYZReader reader(root + "/examples/idr_contact_switching/switching.xyz", 3);
    Frame frame;
    check(reader.next(frame) && frame.coordinates.size() == 3, "XYZ first frame");
    check(close(frame.coordinates[1][0], 0.3), "XYZ angstrom to nm");
    check(reader.next(frame) && close(frame.time_ps, 1000.0), "XYZ explicit time");
#ifdef ENSEMBLEQL_HAS_CHEMFILES
    check(chemfiles_backend_available(), "chemfiles availability flag enabled");
#else
    check(!chemfiles_backend_available(), "chemfiles availability flag disabled");
    check_throws([&] {
        (void)Trajectory::from_files(root + "/tests/data/missing.xtc",
                                     root + "/examples/idr_contact_switching/switching.pdb");
    }, "compressed format reports missing optional backend");
#endif

    const auto model = simple_topology();
    check(Selection("resid 17").resolve(model).size() == 2, "resid selection");
    check(Selection("resid 17:42").resolve(model).size() == 3, "resid range");
    check(Selection("name CA").resolve(model).size() == 3, "name selection");
    check(Selection("resname ARG").resolve(model).size() == 2, "resname selection");
    check(Selection("chain B").resolve(model).size() == 1, "chain selection");
    check(Selection("protein").resolve(model).size() == 4, "protein excludes water");
    check(Selection("resid 17 and name CA").resolve(model).size() == 1, "selection AND");
    check(Selection("resname ARG or resname GLU").resolve(model).size() == 3, "selection OR");
    check(Selection("protein and not hydrogen").resolve(model).size() == 4, "selection NOT");
    check(Selection("name CA or name CB and chain B").resolve(model).size() == 3, "selection AND precedence");
    check(Selection("(name CA or name CB) and chain B").resolve(model).size() == 1, "selection parentheses");
    check_throws([&] { (void)Selection("resid foo").resolve(model); }, "bad resid error");
    check_throws([&] { (void)Selection("resid 999").resolve(model); }, "zero match error");
    check_throws([&] { (void)Selection("resid 17 and").resolve(model); }, "malformed boolean selection");
}

void test_geometry() {
    Frame frame{0.0, {{0.0, 0.0, 0.0}, {0.3, 0.4, 0.0}, {1.0, 0.0, 0.0}}, std::nullopt};
    check(close(distance(frame.coordinates[0], frame.coordinates[1]), 0.5), "distance");
    check(close(minimum_distance(frame, {0}, {1, 2}), 0.5), "minimum distance");
    check(contacts(frame, {0}, {1}, 0.5).size() == 1, "contact cutoff inclusive");
    check(contacts(frame, {0}, {1}, 0.499).empty(), "contact cutoff below boundary");
    Frame rg_frame{0.0, {{-1.0, 0.0, 0.0}, {1.0, 0.0, 0.0}}, std::nullopt};
    check(close(radius_of_gyration(rg_frame, {0, 1}), 1.0), "radius of gyration");
    check(close(minimum_image_distance({0.1, 0.0, 0.0}, {1.9, 0.0, 0.0}, {2.0, 2.0, 2.0}), 0.2),
          "orthorhombic minimum image distance");
    frame.coordinates = {{0.1, 0.0, 0.0}, {1.9, 0.0, 0.0}};
    frame.box_nm = Vec3{2.0, 2.0, 2.0};
    check(close(minimum_distance(frame, {0}, {1}), 0.2), "frame minimum distance uses PBC");
    check(contacts(frame, {0}, {1}, 0.2).size() == 1, "contact crosses periodic boundary");
    check_throws([&] { (void)radius_of_gyration(frame, {0, 1}); }, "periodic RG requires unwrapping");
    check_throws([&] { (void)minimum_image_distance({0, 0, 0}, {1, 0, 0}, {0, 2, 2}); }, "invalid box rejected");
}

void test_periodic_xyz_and_query() {
    const std::string root = ENSEMBLEQL_SOURCE_DIR;
    XYZReader reader(root + "/tests/data/pbc.xyz", 2);
    Frame frame;
    check(reader.next(frame) && frame.box_nm.has_value() && close((*frame.box_nm)[0], 2.0),
          "custom XYZ box metadata");
    check(close(minimum_distance(frame, {0}, {1}), 0.2), "custom XYZ box drives minimum image");
    check(reader.next(frame) && frame.box_nm.has_value() && close((*frame.box_nm)[1], 2.0),
          "extended XYZ Lattice metadata");
    check(close(minimum_distance(frame, {0}, {1}), 0.4), "extended XYZ lattice drives minimum image");

    auto trajectory = Trajectory::from_files(root + "/tests/data/pbc.xyz", root + "/tests/data/pbc.pdb");
    const auto events = Engine().query(trajectory, "FIND CONTACT(resid 1, resid 2, cutoff=0.4nm);");
    check(events.size() == 1 && close(events[0].start_time, 0.0) && close(events[0].end_time, 1.0),
          "streaming periodic contact query");

    XYZReader triclinic(root + "/tests/data/triclinic.xyz", 2);
    check_throws([&] { (void)triclinic.next(frame); }, "triclinic XYZ fails explicitly");
}

void test_events_and_temporal() {
    EventExtractor extractor("CONTACT");
    std::vector<Event> events;
    for (int i = 0; i < 6; ++i) if (auto event = extractor.push(i, i >= 2 && i <= 4)) events.push_back(*event);
    if (auto event = extractor.finish()) events.push_back(*event);
    check(events.size() == 1 && close(events[0].start_time, 2.0) && close(events[0].end_time, 4.0), "contiguous interval");
    check(filter_for(events, 2.0).size() == 1 && filter_for(events, 2.001).empty(), "FOR inclusive boundary");

    const Event a{0, 2, "A", {}, {}}, touching{2, 4, "B", {}, {}}, overlap{1, 3, "B", {}, {}}, late{3.001, 5, "B", {}, {}};
    check(temporal::before(a, touching), "BEFORE includes touching boundary");
    check(temporal::after(touching, a), "AFTER includes touching boundary");
    check(temporal::overlaps(a, touching), "OVERLAPS includes touching boundary");
    check(temporal::overlaps(a, overlap), "OVERLAPS interval");
    check(temporal::followed_by({a}, {touching}, 0.0).size() == 1, "FOLLOWED_BY zero gap");
    check(temporal::followed_by({a}, {late}, 1.0).empty(), "WITHIN rejects over boundary");
    check(temporal::followed_by({a}, {Event{3, 5, "B", {}, {}}}, 1.0).size() == 1, "WITHIN includes boundary");
}

void test_parser_planner_engine() {
    Parser parser;
    auto simple = parser.parse("FIND CONTACT(resid 17, resid 42);");
    check(simple.root->kind == ast::Kind::Contact, "parse simple contact");
    auto duration = parser.parse("FIND CONTACT(resid 17, resid 42) FOR >= 2ns;");
    check(duration.root->kind == ast::Kind::For, "parse FOR");
    auto temporal_query = parser.parse("FIND CONTACT(resid 17, resid 42) FOLLOWED_BY CONTACT(resid 17, resid 53) WITHIN 3ns;");
    check(temporal_query.root->kind == ast::Kind::FollowedBy, "parse FOLLOWED_BY");
    auto nested = parser.parse("FIND (CONTACT(resid 17, resid 42) FOLLOWED_BY CONTACT(resid 17, resid 53) WITHIN 3ns) BEFORE CONTACT(resid 17, resid 42);");
    check(nested.root->kind == ast::Kind::Before, "parse nested temporal expression");
    const auto cutoff_query = parser.parse("FIND CONTACT(resid 17, resid 42, cutoff=0.35nm);");
    check(close(static_cast<const ast::ContactExpr&>(*cutoff_query.root).cutoff.nm, 0.35), "parse contact cutoff");
    const auto precedence = parser.parse("FIND CONTACT(resid 17, resid 42) OR CONTACT(resid 17, resid 53) AND CONTACT(resid 42, resid 53);");
    check(precedence.root->kind == ast::Kind::Or &&
          static_cast<const ast::OrExpr&>(*precedence.root).right->kind == ast::Kind::And,
          "query AND precedence");
    check_throws<QueryError>([&] { (void)parser.parse("FIND CONTACT(resid 17);"); }, "CONTACT arity error");
    check_throws<QueryError>([&] { (void)parser.parse("FIND CONTACT(resid 17, resid 42) FOLLOW CONTACT(resid 1, resid 2);"); }, "unknown operator error");
    check_throws<QueryError>([&] { (void)parser.parse("FIND RG(protein) < 2ns;"); }, "unit dimension error");
    check_throws<QueryError>([&] { (void)parser.parse("FIND CONTACT(resid 17, resid 42, cutoff=2ns);"); }, "contact cutoff dimension error");
    check_throws<QueryError>([&] { (void)parser.parse("FIND CONTACT(resid 17, resid 42) < 1nm;"); }, "boolean contact comparison error");

    const std::string root = ENSEMBLEQL_SOURCE_DIR;
    auto trajectory = Trajectory::from_files(root + "/examples/idr_contact_switching/switching.xyz",
                                             root + "/examples/idr_contact_switching/switching.pdb");
    const auto plan = Planner().plan(temporal_query, trajectory.topology());
    check(plan.unique_observables == 2 && plan.required_selections.size() == 3, "planner requirements and deduplication");
    const auto explanation = Planner().explain(plan);
    check(explanation.streaming && explanation.observables.size() == 2, "plan explanation streaming and observables");
    check(explanation.frame_predicates.size() == 2 && explanation.temporal_operations.size() == 1,
          "plan explanation stages");
    check(explanation.temporal_operations[0] == "FOLLOWED_BY WITHIN 3000ps", "plan explanation temporal detail");
    check(explanation.tree.find("FOLLOWED_BY") != std::string::npos, "plan explanation tree");
    const auto results = Engine().execute(trajectory, plan);
    check(results.size() == 1, "switching event found");
    if (!results.empty()) {
        check(close(results[0].start_time, 0.0) && close(results[0].end_time, 4000.0), "switching event boundaries");
        check(results[0].metadata.at("transition_gap_ps") == "0.000000", "switching transition gap");
    }
    const auto persistent = Engine().query(trajectory, "FIND CONTACT(resid 17, resid 42) FOR >= 2ns;");
    check(persistent.size() == 1 && close(persistent[0].duration(), 2000.0), "streaming query with FOR");
    const auto count_events = Engine().query(trajectory, "FIND CONTACT_COUNT(resid 17, resid 42) >= 1;");
    check(count_events.size() == 1 && close(count_events[0].end_time, 2000.0), "CONTACT_COUNT comparison execution");
    const auto tighter = Engine().query(trajectory, "FIND CONTACT(resid 17, resid 42, cutoff=0.35nm);");
    check(tighter.size() == 1 && close(tighter[0].end_time, 1000.0), "configurable inclusive contact cutoff");
    const auto shared = parser.parse("FIND DISTANCE(resid 17, resid 42) < 1nm AND DISTANCE(resid 42, resid 17) < 2nm;");
    check(Planner().plan(shared, trajectory.topology()).unique_observables == 1, "planner shares symmetric observable");
    check_throws<QueryError>([&] { (void)Planner().plan(parser.parse("FIND RG(protein);"), trajectory.topology()); }, "numeric root requires comparison");
    check(Engine().explain(trajectory.topology(), "FIND CONTACT(resid 17, resid 42);").frame_predicates.size() == 1,
          "engine explain without trajectory scan");
}

void test_sampled_time_semantics() {
    Topology topology({{0, "CA", "ARG", 1, 'A', "C"}, {1, "CA", "ASP", 2, 'A', "C"}});
    const auto frame = [](double time, double separation) {
        return Frame{time, {{0.0, 0.0, 0.0}, {separation, 0.0, 0.0}}, std::nullopt};
    };
    auto irregular_reader = std::make_shared<VectorReader>(std::vector<Frame>{frame(0, 0.3), frame(3, 0.3), frame(10, 0.3)});
    Trajectory irregular(topology, irregular_reader);
    const auto spanning = Engine().query(irregular, "FIND CONTACT(resid 1, resid 2) FOR >= 10ps;");
    check(spanning.size() == 1 && close(spanning[0].duration(), 10.0), "irregular timestamps use observed span");

    auto one_reader = std::make_shared<VectorReader>(std::vector<Frame>{frame(0, 1.0), frame(3, 0.3), frame(10, 1.0)});
    Trajectory one_sample(topology, one_reader);
    const auto instantaneous = Engine().query(one_sample, "FIND CONTACT(resid 1, resid 2);");
    check(instantaneous.size() == 1 && close(instantaneous[0].duration(), 0.0), "single-sample event has zero observed duration");
    check(Engine().query(one_sample, "FIND CONTACT(resid 1, resid 2) FOR >= 0.001ps;").empty(), "FOR rejects single-sample event above zero");

    auto invalid_reader = std::make_shared<VectorReader>(std::vector<Frame>{frame(1, 0.3), frame(1, 0.3)});
    Trajectory invalid(topology, invalid_reader);
    check_throws([&] { (void)Engine().query(invalid, "FIND CONTACT(resid 1, resid 2);"); }, "non-monotonic timestamps rejected");
}
} // namespace

int main() {
    test_units();
    test_io_and_selections();
    test_geometry();
    test_periodic_xyz_and_query();
    test_events_and_temporal();
    test_parser_planner_engine();
    test_sampled_time_semantics();
    if (failures) { std::cerr << failures << " test(s) failed\n"; return 1; }
    std::cout << "All EnsembleQL C++ tests passed\n";
    return 0;
}
