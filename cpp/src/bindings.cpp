#include "ensembleql/engine.hpp"
#include "ensembleql/geometry.hpp"
#include "ensembleql/parser.hpp"
#include "ensembleql/selection.hpp"
#include "ensembleql/topology.hpp"
#include "ensembleql/trajectory.hpp"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;
using namespace ensembleql;

namespace {
py::dict explanation_dict(const PlanExplanation& explanation) {
    py::dict result;
    result["streaming"] = explanation.streaming;
    result["selections"] = explanation.selections;
    result["observables"] = explanation.observables;
    result["frame_predicates"] = explanation.frame_predicates;
    result["temporal_operations"] = explanation.temporal_operations;
    result["plan_tree"] = explanation.tree;
    return result;
}
} // namespace

PYBIND11_MODULE(_core, module) {
    module.doc() = "EnsembleQL C++20 execution engine";
    py::register_exception<QueryError>(module, "QueryError");

    py::class_<Atom>(module, "Atom")
        .def_readonly("index", &Atom::index).def_readonly("name", &Atom::name)
        .def_readonly("residue_name", &Atom::residue_name).def_readonly("residue_index", &Atom::residue_index)
        .def_readonly("chain", &Atom::chain).def_readonly("element", &Atom::element);
    py::class_<Topology>(module, "Topology")
        .def_static("from_pdb", &Topology::from_pdb).def_property_readonly("atoms", &Topology::atoms)
        .def("select", [](const Topology& topology, const std::string& expression) { return Selection(expression).resolve(topology); });
    py::class_<Event>(module, "Event")
        .def_readonly("start", &Event::start_time).def_readonly("end", &Event::end_time)
        .def_property_readonly("duration", &Event::duration).def_readonly("type", &Event::type)
        .def_readonly("selections", &Event::selections).def_readonly("metadata", &Event::metadata)
        .def("__repr__", [](const Event& event) { return "Event(type='" + event.type + "', start=" + std::to_string(event.start_time) + ", end=" + std::to_string(event.end_time) + ")"; });
    py::class_<Trajectory>(module, "NativeTrajectory")
        .def_static("from_files", &Trajectory::from_files)
        .def_property_readonly("topology", &Trajectory::topology, py::return_value_policy::reference_internal)
        .def("query", [](Trajectory& trajectory, const std::string& text) { return Engine().query(trajectory, text); })
        .def("explain", [](const Trajectory& trajectory, const std::string& text) {
            return explanation_dict(Engine().explain(trajectory.topology(), text));
        });

    module.def("parse_query", [](const std::string& text) {
        const auto query = Parser().parse(text);
        return static_cast<int>(query.root->kind);
    });
    module.def("distance", [](const Vec3& a, const Vec3& b) { return distance(a, b); });
    module.def("minimum_image_distance", &minimum_image_distance,
               py::arg("a"), py::arg("b"), py::arg("box_nm"));
    module.def("chemfiles_backend_available", &chemfiles_backend_available);
    module.def("explain_query", [](const Topology& topology, const std::string& text) {
        return explanation_dict(Engine().explain(topology, text));
    });
}
