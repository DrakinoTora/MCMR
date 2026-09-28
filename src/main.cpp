#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "simulation.hpp"
#include "simulation_mg.hpp"
#include "material.hpp"

namespace py = pybind11;

PYBIND11_MODULE(_mcmr_cpp, m) {
    // Single source of truth for material name resolution: the SAME alias
    // table (via get_material_info) that Simulation/SimulationMG use
    // internally to build the geometry. Exposed so Python-side code (group
    // materials, validation) can normalize names to the exact canonical
    // spelling the C++ engines key their per-material data on, instead of
    // duplicating -- and risking drifting from -- the alias table.
    m.def("canonical_material_name", [](const std::string& name) {
        return get_material_info(name).symbol;
    }, py::arg("name"),
       "Resolve a material name/alias (e.g. 'fe', 'besi', 'iron') to its "
       "canonical symbol ('Fe'). Raises ValueError if the name is unknown.");

    py::class_<Tally>(m, "Tally")
        .def_readonly("transmission", &Tally::transmission)
        .def_readonly("time_taken", &Tally::time_taken)
        .def_readonly("E_born", &Tally::E_born)
        .def_readonly("E_leak", &Tally::E_leak)
        .def_readonly("G_born", &Tally::G_born)
        .def_readonly("G_leak", &Tally::G_leak)
        .def_readonly("x_history", &Tally::x_history)
        .def_readonly("y_history", &Tally::y_history)
        .def_readonly("history_generation", &Tally::history_generation)
        .def_readonly("k_used", &Tally::k_used)
        .def_readonly("k_estimate", &Tally::k_estimate)
        // Region-based tally grid (replaces the old absorp_by_material /
        // fission_by_material per-material counters): a resolution grid laid
        // independently over the whole world -- every fission/absorption event
        // adds 1 point to the cell it occurred in. region_tally is flat,
        // row-major (index = iy_t*tally_nx + ix_t, iy_t=0 = bottom row/y=0,
        // ix_t=0 = leftmost col/x=0); reshape client-side using tally_nx/tally_ny.
        .def_readonly("fission_tally", &Tally::fission_tally)
        .def_readonly("absorp_tally", &Tally::absorp_tally)
        .def_readonly("tally_nx", &Tally::tally_nx)
        .def_readonly("tally_ny", &Tally::tally_ny)
        .def_readonly("tally_dx", &Tally::tally_dx)
        .def_readonly("tally_dy", &Tally::tally_dy);

    py::class_<Simulation>(m, "Simulation")
        .def(py::init<int, double, double,
                    const std::vector<double>&,
                    const std::vector<double>&,
                    const std::vector<std::vector<std::string>>&,
                    const std::vector<std::vector<double>>&,
                    const std::vector<double>&,
                    const std::vector<double>&,
                    const std::vector<double>&,
                    const std::vector<std::string>&,
                    const std::vector<double>&,
                    int,
                    int,
                    int,
                    const std::string&,
                    const std::string&,
                    const std::string&,
                    const std::string&>(),
            py::arg("N"),
            py::arg("x_world"), py::arg("y_world"),
            py::arg("x_grid"), py::arg("y_grid"),
            py::arg("material_matrix"),
            py::arg("sources"),
            py::arg("circle_cx") = std::vector<double>{},
            py::arg("circle_cy") = std::vector<double>{},
            py::arg("circle_r") = std::vector<double>{},
            py::arg("circle_material") = std::vector<std::string>{},
            py::arg("circle_source") = std::vector<double>{},
            py::arg("max_history_save") = 50,
            py::arg("tally_res_x") = 0,
            py::arg("tally_res_y") = 0,
            py::arg("bc_top") = "vacuum",
            py::arg("bc_bot") = "vacuum",
            py::arg("bc_left") = "vacuum",
            py::arg("bc_right") = "vacuum")
        .def("set_cross_sections", &Simulation::set_cross_sections)
        .def("run", &Simulation::run)
        .def("export_xml", &Simulation::export_xml, py::arg("filename") = "mcmr_results.xml")
        .def("get_tally", &Simulation::get_tally);

    py::class_<SimulationMG>(m, "SimulationMG")
        .def(py::init<int, double, double,
                    const std::vector<double>&,
                    const std::vector<double>&,
                    const std::vector<std::vector<std::string>>&,
                    const std::vector<std::vector<double>>&,
                    const std::vector<double>&,
                    const std::vector<double>&,
                    const std::vector<double>&,
                    const std::vector<std::string>&,
                    const std::vector<double>&,
                    int,
                    int,
                    int,
                    const std::string&,
                    const std::string&,
                    const std::string&,
                    const std::string&>(),
            py::arg("N"),
            py::arg("x_world"), py::arg("y_world"),
            py::arg("x_grid"), py::arg("y_grid"),
            py::arg("material_matrix"),
            py::arg("sources"),
            py::arg("circle_cx") = std::vector<double>{},
            py::arg("circle_cy") = std::vector<double>{},
            py::arg("circle_r") = std::vector<double>{},
            py::arg("circle_material") = std::vector<std::string>{},
            py::arg("circle_source") = std::vector<double>{},
            py::arg("max_history_save") = 10,
            py::arg("tally_res_x") = 0,
            py::arg("tally_res_y") = 0,
            py::arg("bc_top") = "vacuum",
            py::arg("bc_bot") = "vacuum",
            py::arg("bc_left") = "vacuum",
            py::arg("bc_right") = "vacuum")
        .def("set_group_data", &SimulationMG::set_group_data,
            py::arg("n_groups"), py::arg("sigma_t"), py::arg("sigma_s"),
            py::arg("sigma_f"), py::arg("nu"))
        .def("run", &SimulationMG::run, py::arg("n_gen") = 1, py::arg("k") = 1.0)
        .def("export_xml", &SimulationMG::export_xml, py::arg("filename") = "mcmr_results_mg.xml")
        .def("get_tally", &SimulationMG::get_tally);
}