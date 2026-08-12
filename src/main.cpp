#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "simulation.hpp"

namespace py = pybind11;

PYBIND11_MODULE(_mcmr_cpp, m) {
    py::class_<Tally>(m, "Tally")
        .def_readonly("absorp_by_material", &Tally::absorp_by_material)
        .def_readonly("transmisi", &Tally::transmisi)
        .def_readonly("time_taken", &Tally::time_taken)
        .def_readonly("E_born", &Tally::E_born)
        .def_readonly("E_leak", &Tally::E_leak)
        .def_readonly("x_history", &Tally::x_history)
        .def_readonly("y_history", &Tally::y_history);

    py::class_<Simulation>(m, "Simulation")
        .def(py::init<int, double, double,
                       const std::vector<double>&,
                       const std::vector<double>&,
                       const std::vector<std::vector<std::string>>&,
                       int>(),
             py::arg("N"),
             py::arg("x_world"), py::arg("y_world"),
             py::arg("x_grid"), py::arg("y_grid"),
             py::arg("material_matrix"),
             py::arg("max_history_save") = 50)
        .def("set_cross_sections", &Simulation::set_cross_sections)
        .def("run", &Simulation::run)
        .def("export_xml", &Simulation::export_xml, py::arg("filename") = "mcmr_results.xml")
        .def("get_tally", &Simulation::get_tally);
}
