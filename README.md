# MCMR (Monte Carlo Multi-Region)
MCMR is a high-performance C++17 Monte Carlo Neutron Transport simulation library, wrapped for Python via pybind11 and scikit-build-core. It simulates 2D neutron transport, scattering, absorption, fission, and transmission across multi-region geometries built from any combination of `Be`, `C`, `Fe`, and `Pb` regions, using analog techniques.

## Features
- **High Performance**: core Monte Carlo tracking engine written in C++17.
- **Pythonic API**: `World` and a `Geometry` builder for defining regions, materials, and sources without hand-building grids/matrices; both produce the same `World` object, which drives the simulation via pybind11.
- **Rectangular + Circular Regions**: build a world out of a rectangular grid, then layer true circular regions on top with `World.add_circle()` / `World.add_circles()` (or the matching methods on `Geometry`) — circles use real line-circle intersection in the physics engine, not a rasterized approximation.
- **Two Physics Modes**:
  - `mode="continuous"` (default) — Watt fission spectrum birth + continuous-energy cross sections, interpolated from the bundled nuclear data.
  - `mode="group"` — multi-group (discretized energy) transport, including fission ($\Sigma_f$, $\nu$), driven by user-supplied `GroupMaterial` / `MaterialLibrary` cross-section sets.
- **Flexible Boundary Conditions**: `vacuum` or `reflective`, set independently per side (top/bottom/left/right).
- **Self-Contained Cross-Section Data**: built-in nuclear cross-section HTML data parser (BeautifulSoup4 + pandas), with fast linear interpolation for macroscopic cross-sections ($\Sigma_t$ and $\Sigma_s$).
- **XML Import/Export**: save a `World` definition, a `MaterialLibrary` (group cross sections), and simulation results (particle tallies, trajectories) all to XML, and reload any of them later to re-run or re-plot.
- **Trajectory Visualization**: `ResultsPlotter` renders neutron trajectories over the region geometry (rectangles and circles alike) with matplotlib.

## Prerequisites
Before building and installing the library, ensure you have the following installed on your system:
- C++ Compiler: GCC/Clang with C++17 support.
- CMake: Version 3.15 or higher.
- Python: Version 3.8 or higher.
- Python Packages: pip, wheel.

On Ubuntu/Debian, you can install the build dependencies using:

```bash
sudo apt update
sudo apt install build-essential cmake python3-dev
```

## Installation
1. Clone the Repository

```bash
git clone https://github.com/DrakinoTora/MCMR.git
cd MCMR
```

2. Install Required Python Dependencies

```bash
pip install scikit-build-core pybind11 pandas beautifulsoup4 matplotlib numpy
```

3. Build and Install MCMR

```bash
pip install .
```

## Testing the Installation
You can quickly verify that the library is installed properly by running the following command in your terminal:

```bash
python3 -c "import mcmr; print(mcmr.__all__)"
```

Expected Output:

```plaintext
['Tally', 'Simulation', 'SimulationMG', 'ResultsPlotter', 'World', 'Geometry', 'GroupMaterial', 'MaterialLibrary', 'load_all_materials', 'load_cross_section']
```
