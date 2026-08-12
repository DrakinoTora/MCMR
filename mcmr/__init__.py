from ._mcmr_cpp import Tally, Simulation
from .cross_section import load_cross_section, load_all_materials
from .plotter import ResultsPlotter


def load_and_run(N, x_world, y_world, x_grid, y_grid, material_matrix, max_save=50):
    """Convenience wrapper ala OpenMC.

    x_world, y_world : ukuran dunia simulasi (0,0) sampai (x_world, y_world)
    x_grid, y_grid    : batas interior grid, terurut naik, di antara 0 dan x_world/y_world
                         mis. x_grid=[3,15,21] -> 4 sel arah x: [0,3],[3,15],[15,21],[21,x_world]
    material_matrix   : list 2D nama material, ukurannya harus
                         (len(x_grid)+1) baris x (len(y_grid)+1) kolom
                         material_matrix[i][j] = material di sel kolom-x ke-i, baris-y ke-j
    """
    sim = Simulation(N, x_world, y_world, x_grid, y_grid, material_matrix, max_save)
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
