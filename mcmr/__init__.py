from ._mcmr_cpp import Tally, Simulation
from .cross_section import load_cross_section, load_all_materials
from .plotter import ResultsPlotter
from .world import export_world, load_world


def load_and_run(N, x_world, y_world, x_grid, y_grid, material_matrix, sources, max_save=50,
                  bc_top="vacuum", bc_bot="vacuum", bc_left="vacuum", bc_right="vacuum"):
    """Convenience wrapper ala OpenMC.

    ...
    material_matrix, sources : list 2D berukuran (len(y_grid)+1) baris x (len(x_grid)+1) kolom.
                         Konvensi index [row][col]: row=0 adalah baris PALING ATAS (y tertinggi),
                         menurun sampai row terakhir (y=0). col=0 adalah kolom PALING KIRI (x=0),
                         naik sampai col terakhir (x tertinggi). Yaitu ditulis natural seperti
                         menggambar grid di kertas. `sources` punya ukuran & orientasi yang sama.
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


def run_from_world(N, world_xml, max_save=50):
    """Jalankan simulasi langsung dari file world XML hasil mcmr.export_world().

    Menghindari kerja dua kali: definisikan dunia sekali lewat export_world(),
    lalu cukup panggil run_from_world(N, "world.xml") -- tidak perlu ketik ulang
    x_world, y_world, x_grid, y_grid, material_matrix, sources, bc_*.

    N         : jumlah partikel neutron yang disimulasikan
    world_xml : path file world XML (hasil export_world)
    max_save  : jumlah maksimum lintasan neutron yang disimpan untuk plotting
    """
    w = load_world(world_xml)
    return load_and_run(
        N=N,
        x_world=w["x_world"], y_world=w["y_world"],
        x_grid=w["x_grid"], y_grid=w["y_grid"],
        material_matrix=w["material_matrix"], sources=w["sources"],
        max_save=max_save,
        bc_top=w["bc_top"], bc_bot=w["bc_bot"],
        bc_left=w["bc_left"], bc_right=w["bc_right"],
    )


__all__ = [
    "Tally",
    "Simulation",
    "ResultsPlotter",
    "load_and_run",
    "run_from_world",
    "export_world",
    "load_world",
    "load_all_materials",
]