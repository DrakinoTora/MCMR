#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "simulation.hpp"

namespace py = pybind11;

PYBIND11_MODULE(_mcmr_cpp, m) {
    py::class_<BoxBoundary>(m, "BoxBoundary")
        .def(py::init<double, double, double, double>(), py::arg("x1"), py::arg("y1"), py::arg("x2"), py::arg("y2"));

    py::class_<WorldBoundary>(m, "WorldBoundary")
        .def(py::init<double, double>(), py::arg("max_x"), py::arg("max_y"));

    py::class_<Tally>(m, "Tally")
        .def_readonly("absorp_material", &Tally::absorp_material)
        .def_readonly("absorp_outside", &Tally::absorp_outside)
        .def_readonly("transmisi", &Tally::transmisi)
        .def_readonly("time_taken", &Tally::time_taken)
        .def_readonly("E_born", &Tally::E_born)
        .def_readonly("E_leak", &Tally::E_leak)
        .def_readonly("x_history", &Tally::x_history)
        .def_readonly("y_history", &Tally::y_history);

    py::class_<Simulation>(m, "Simulation")
        .def(py::init<int, const std::string&, BoxBoundary, WorldBoundary, int>(),
             py::arg("N"), py::arg("material_name"), py::arg("box"), py::arg("world"), py::arg("max_history_save") = 50)
        .def("set_cross_sections", &Simulation::set_cross_sections)
        .def("run", &Simulation::run)
        .def("export_xml", &Simulation::export_xml, py::arg("filename") = "mcmr_results.xml")
        .def("get_tally", &Simulation::get_tally);
}