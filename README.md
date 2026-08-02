# MCMR (Monte Carlo Multi-Region)
MCMR is a high-performance C++ Monte Carlo Neutron Transport simulation library wrapped for Python using pybind11 and scikit-build-core. It simulates 2D neutron transport, scattering, absorption, and transmission across multi-region geometries ($Fe - Pb - Fe$) using both analog and non-analog (implicit capture with Russian Roulette) techniques.

## Features
- High Performance: Core Monte Carlo tracking engine written in C++17.
- Seamless Python Binding: Exposed directly to Python via pybind11.
- Self-Contained Data: Internal nuclear cross-section HTML data parser (BeautifulSoup4 + pandas).
- Cross-Section Interpolation: Fast linear interpolation for macroscopic cross-sections ($\Sigma_t$ and $\Sigma_s$).
- Fission Energy Spectrum: Built-in Watt Fission Spectrum sampling.

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
git clone https://github.com/your-username/MCMR.git
cd MCMR
```

2. Install Required Python Dependencies

```bash
pip install scikit-build-core pybind11 pandas beautifulsoup4
```

3. Build and Install MCMR

```bash
pip install .
```

## Usage in Python
Once installed, you can import and run the mcmr library directly inside any Python script or Jupyter Notebook (.ipynb).

## Testing the Installation
You can quickly verify that the library is installed properly by running the following command in your terminal:

```bash
python3 -c "import mcmr; print(dir(mcmr))"
```

Expected Output:

```plaintext
['SimulationResult', '__all__', 'load_cross_section', 'neutron_sim', 'simpan_hasil_xml']
```