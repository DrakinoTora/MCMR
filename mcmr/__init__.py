from ._mcmr_cpp import Tally, Simulation
from .cross_section import load_cross_section, load_all_materials
from .plotter import ResultsPlotter


def load_and_run(N, x_world, y_world, x_grid, y_grid, material_matrix, sources, max_save=50,
                  bc_top="vacuum", bc_bot="vacuum", bc_left="vacuum", bc_right="vacuum"):
    """Convenience wrapper ala OpenMC.

    ...
    bc_top, bc_bot, bc_left, bc_right : kondisi batas di tepi dunia simulasi.
                         Pilihan: "vacuum" (neutron mati/leak, default) atau
                         "reflective" (neutron memantul balik, energi tetap).
                         Sisi yang tidak didefinisikan otomatis "vacuum".
    """
    sim = Simulation(N, x_world, y_world, x_grid, y_grid, material_matrix, sources, max_save,
                      bc_top, bc_bot, bc_left, bc_right)
    E_tot, Sig_tot, E_scat, Sig_scat = load_all_materials()
    sim.set_cross_sections(E_tot, Sig_tot, E_scat, Sig_scat)

    sim.run()
    return sim


__all__ = [
    "Tally",
    "Simulation",
    "ResultsPlotter",
    "load_and_run",
    "load_all_materials",
]
