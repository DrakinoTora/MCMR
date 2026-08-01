#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "simulation.hpp"
#include "physics.hpp"
#include "io.hpp"

namespace py = pybind11;

PYBIND11_MODULE(_mcmr_cpp, m) {
    m.doc() = "Library Monte Carlo Transport Neutron (MCMR) C++ Python Binding";

    // Bind Struct Result
    py::class_<SimulationResult>(m, "SimulationResult")
        .def_readonly("abs_A", &SimulationResult::abs_A)
        .def_readonly("abs_B", &SimulationResult::abs_B)
        .def_readonly("transmisi", &SimulationResult::transmisi)
        .def_readonly("time_taken", &SimulationResult::time_taken)
        .def_readonly("E_start", &SimulationResult::E_start)
        .def_readonly("E_A_escape", &SimulationResult::E_A_escape)
        .def_readonly("E_B_escape", &SimulationResult::E_B_escape);

    // Bind Fungsi Simulasi Utama
    m.def("neutron_sim", &neutron_sim, "Jalankan Simulasi Monte Carlo Neutron",
          py::arg("a"), py::arg("b"), py::arg("c"), py::arg("N_partikel"), py::arg("analog"),
          py::arg("E_Al_t"), py::arg("Sig_Al_t"), py::arg("E_Al_s"), py::arg("Sig_Al_s"),
          py::arg("E_Pb_t"), py::arg("Sig_Pb_t"), py::arg("E_Pb_s"), py::arg("Sig_Pb_s"));

    // Bind Fungsi Pembantu
    m.def("simpan_hasil_xml", &simpan_hasil_xml, "Simpan hasil ke file XML",
          py::arg("res"), py::arg("filename") = "hasil_simulasi.xml");
}