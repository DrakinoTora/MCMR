#include "material.hpp"
#include <stdexcept>
#include <algorithm>
#include <unordered_map>

namespace {

// Full periodic table (Z = 1..118): canonical symbol, English name, atomic mass.
// Mass is the standard atomic weight; for elements with no stable isotope it is
// the mass number of the longest-lived isotope. It is only informational here --
// elements outside Be/C/Fe/Pb have no continuous-energy data, so this value is
// not used for transport (multi-group mode keys everything on the symbol).
struct Element {
    const char* symbol;
    const char* name;
    double A;
};

const Element ELEMENTS[] = {
    {"H",  "hydrogen",      1.008},        {"He", "helium",        4.002602},
    {"Li", "lithium",       6.94},         {"Be", "beryllium",     9.0121831},
    {"B",  "boron",         10.81},        {"C",  "carbon",        12.011},
    {"N",  "nitrogen",      14.007},       {"O",  "oxygen",        15.999},
    {"F",  "fluorine",      18.998403163}, {"Ne", "neon",          20.1797},
    {"Na", "sodium",        22.98976928},  {"Mg", "magnesium",     24.305},
    {"Al", "aluminium",     26.9815385},   {"Si", "silicon",       28.085},
    {"P",  "phosphorus",    30.973761998}, {"S",  "sulfur",        32.06},
    {"Cl", "chlorine",      35.45},        {"Ar", "argon",         39.948},
    {"K",  "potassium",     39.0983},      {"Ca", "calcium",       40.078},
    {"Sc", "scandium",      44.955908},    {"Ti", "titanium",      47.867},
    {"V",  "vanadium",      50.9415},      {"Cr", "chromium",      51.9961},
    {"Mn", "manganese",     54.938044},    {"Fe", "iron",          55.845},
    {"Co", "cobalt",        58.933194},    {"Ni", "nickel",        58.6934},
    {"Cu", "copper",        63.546},       {"Zn", "zinc",          65.38},
    {"Ga", "gallium",       69.723},       {"Ge", "germanium",     72.630},
    {"As", "arsenic",       74.921595},    {"Se", "selenium",      78.971},
    {"Br", "bromine",       79.904},       {"Kr", "krypton",       83.798},
    {"Rb", "rubidium",      85.4678},      {"Sr", "strontium",     87.62},
    {"Y",  "yttrium",       88.90584},     {"Zr", "zirconium",     91.224},
    {"Nb", "niobium",       92.90637},     {"Mo", "molybdenum",    95.95},
    {"Tc", "technetium",    98.0},         {"Ru", "ruthenium",     101.07},
    {"Rh", "rhodium",       102.90550},    {"Pd", "palladium",     106.42},
    {"Ag", "silver",        107.8682},     {"Cd", "cadmium",       112.414},
    {"In", "indium",        114.818},      {"Sn", "tin",           118.710},
    {"Sb", "antimony",      121.760},      {"Te", "tellurium",     127.60},
    {"I",  "iodine",        126.90447},    {"Xe", "xenon",         131.293},
    {"Cs", "caesium",       132.90545196}, {"Ba", "barium",        137.327},
    {"La", "lanthanum",     138.90547},    {"Ce", "cerium",        140.116},
    {"Pr", "praseodymium",  140.90766},    {"Nd", "neodymium",     144.242},
    {"Pm", "promethium",    145.0},        {"Sm", "samarium",      150.36},
    {"Eu", "europium",      151.964},      {"Gd", "gadolinium",    157.25},
    {"Tb", "terbium",       158.92535},    {"Dy", "dysprosium",    162.500},
    {"Ho", "holmium",       164.93033},    {"Er", "erbium",        167.259},
    {"Tm", "thulium",       168.93422},    {"Yb", "ytterbium",     173.045},
    {"Lu", "lutetium",      174.9668},     {"Hf", "hafnium",       178.49},
    {"Ta", "tantalum",      180.94788},    {"W",  "tungsten",      183.84},
    {"Re", "rhenium",       186.207},      {"Os", "osmium",        190.23},
    {"Ir", "iridium",       192.217},      {"Pt", "platinum",      195.084},
    {"Au", "gold",          196.966569},   {"Hg", "mercury",       200.592},
    {"Tl", "thallium",      204.38},       {"Pb", "lead",          207.2},
    {"Bi", "bismuth",       208.98040},    {"Po", "polonium",      209.0},
    {"At", "astatine",      210.0},        {"Rn", "radon",         222.0},
    {"Fr", "francium",      223.0},        {"Ra", "radium",        226.0},
    {"Ac", "actinium",      227.0},        {"Th", "thorium",       232.0377},
    {"Pa", "protactinium",  231.03588},    {"U",  "uranium",       238.02891},
    {"Np", "neptunium",     237.0},        {"Pu", "plutonium",     244.0},
    {"Am", "americium",     243.0},        {"Cm", "curium",        247.0},
    {"Bk", "berkelium",     247.0},        {"Cf", "californium",   251.0},
    {"Es", "einsteinium",   252.0},        {"Fm", "fermium",       257.0},
    {"Md", "mendelevium",   258.0},        {"No", "nobelium",      259.0},
    {"Lr", "lawrencium",    266.0},        {"Rf", "rutherfordium", 267.0},
    {"Db", "dubnium",       268.0},        {"Sg", "seaborgium",    269.0},
    {"Bh", "bohrium",       270.0},        {"Hs", "hassium",       277.0},
    {"Mt", "meitnerium",    278.0},        {"Ds", "darmstadtium",  281.0},
    {"Rg", "roentgenium",   282.0},        {"Cn", "copernicium",   285.0},
    {"Nh", "nihonium",      286.0},        {"Fl", "flerovium",     289.0},
    {"Mc", "moscovium",     290.0},        {"Lv", "livermorium",   293.0},
    {"Ts", "tennessine",    294.0},        {"Og", "oganesson",     294.0},
};
constexpr std::size_t N_ELEMENTS = sizeof(ELEMENTS) / sizeof(ELEMENTS[0]);
static_assert(N_ELEMENTS == 118, "periodic table must list exactly 118 elements");

std::string to_lower(std::string s) {
    std::transform(s.begin(), s.end(), s.begin(), ::tolower);
    return s;
}

// lowercase symbol / name / alias -> index into ELEMENTS
const std::unordered_map<std::string, std::size_t>& element_index() {
    static const std::unordered_map<std::string, std::size_t> index = [] {
        std::unordered_map<std::string, std::size_t> m;
        for (std::size_t i = 0; i < N_ELEMENTS; ++i) {
            m.emplace(to_lower(ELEMENTS[i].symbol), i);
            m.emplace(to_lower(ELEMENTS[i].name), i);
        }
        // alternative spellings of the English names
        auto alias = [&](const char* spelling, const char* symbol) {
            m.emplace(spelling, m.at(to_lower(symbol)));
        };
        alias("aluminum", "Al");
        alias("sulphur", "S");
        alias("cesium", "Cs");
        // Indonesian names that were already accepted before the full table existed
        alias("berilium", "Be");
        alias("grafit", "C");
        alias("besi", "Fe");
        alias("timbal", "Pb");
        return m;
    }();
    return index;
}

}  // namespace

MaterialInfo get_material_info(const std::string& name) {
    const auto& index = element_index();
    auto it = index.find(to_lower(name));
    if (it == index.end())
        throw std::invalid_argument(
            "Material '" + name + "' unknown. Use any element symbol or English name "
            "(e.g. U, Na, Fe, hydrogen, uranium).");

    const Element& e = ELEMENTS[it->second];
    const std::string sym = e.symbol;

    // The four materials that ship with continuous-energy cross-section data
    // (data/*.html). Atomic number density N = rho * N_A / M [atoms/cm^3], from
    // standard density (g/cm^3) and molar mass (g/mol). N_A = 6.02214076e23 /mol.
    //   Be: rho=1.848,  M=9.012182  -> N=1.234875e23
    //   C  (graphite): rho=2.260, M=12.0107   -> N=1.133159e23
    //   Fe: rho=7.874,  M=55.845   -> N=8.491062e22
    //   Pb: rho=11.34,  M=207.2    -> N=3.295901e22
    if (sym == "Be") return {"Be", 425, 9.0, 1.234875e23};
    if (sym == "C")  return {"C", 600, 12.0, 1.133159e23};
    if (sym == "Fe") return {"Fe", 2631, 56.0, 8.491062e22};
    if (sym == "Pb") return {"Pb", 8237, 208.0, 3.295901e22};

    // Any other element: valid for multi-group mode, where the user supplies
    // the cross sections. mat_code = -1 marks "no continuous-energy data", and
    // Simulation::run() rejects such a material up front with a clear message.
    return {sym, -1, e.A, 0.0};
}
