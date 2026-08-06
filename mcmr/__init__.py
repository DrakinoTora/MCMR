from ._mcmr_cpp import BoxBoundary, WorldBoundary, Tally, Simulation
from .cross_section import load_cross_section, load_all_materials
from .plotter import ResultsPlotter

def load_and_run(N, material_name, box_coords, world_coords, max_save=50):
    """Convenience Wrapper ala OpenMC"""
    box = BoxBoundary(*box_coords)
    world = WorldBoundary(*world_coords)
    
    sim = Simulation(N, material_name, box, world, max_save)
    E_tot, Sig_tot, E_scat, Sig_scat = load_all_materials()
    sim.set_cross_sections(E_tot, Sig_tot, E_scat, Sig_scat)
    
    sim.run()
    return sim

__all__ = [
    "BoxBoundary",
    "WorldBoundary",
    "Tally",
    "Simulation",
    "ResultsPlotter",
    "load_and_run",
    "load_all_materials"
]