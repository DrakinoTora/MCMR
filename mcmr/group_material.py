import math
import xml.etree.ElementTree as ET
from ._mcmr_cpp import canonical_material_name

__all__ = ["GroupMaterial", "MaterialLibrary"]


class GroupMaterial:
    """Multi-group (discretized energy) cross-section data for ONE material.

    `name` can be ANY chemical element, by symbol or English name, in any
    letter case -- "U", "u", "uranium", "Fe", "besi", "iron", etc. It is
    resolved to the element's canonical symbol (e.g. "U", "Fe") right here,
    through the SAME lookup the C++ engine uses internally
    (see get_material_info() in src/material.cpp, exposed to Python as
    canonical_material_name()). This is the single source of truth for
    aliases: nothing else in this file, or in World._run_group(), keeps its
    own copy of the alias table, so a GroupMaterial and the material_matrix
    cell it's meant for always line up in the end, no matter which spelling
    either one used.

    Number of groups (n_groups) is inferred from len(sigma_t) -- there is no
    separate "how many groups" parameter to keep in sync by hand.

        Fe = GroupMaterial(
            name="Fe",
            sigma_t=[0.5, 0.2, 0.1],
            sigma_s=[0.3, 0.1, 0.06],
            sigma_f=[0.001, 0.0005, 0.00001],   # optional, default 0 (no fission)
            nu=[2.7, 2.5, 2.3],                 # optional, default 0
            chi_fission=[0.8, 0.2, 0.0],        # optional, see below
            chi_scatter=[[0.6, 0.3, 0.1],       # optional, see below
                         [0.0, 0.7, 0.3],
                         [0.0, 0.0, 1.0]],
        )

    chi_fission is the fission spectrum: the probability that a neutron born
    from a fission in this material lands in each group. It has one entry per
    group and does NOT depend on the group of the neutron that caused the
    fission. It is normalized internally, so the entries only need to be
    non-negative with a positive sum -- [4, 1, 0] means exactly the same as
    [0.8, 0.2, 0.0]. Every child of a fission draws its own group from it.
    Left out (None or []), fission children inherit the group of the neutron
    that fissioned.

    chi_scatter is the scattering matrix: an n_groups x n_groups list of lists
    (or 2-D array) whose entry [i][j] is the probability that a neutron which
    scatters while in group i continues in group j. So the ROW is the group
    before the scatter and the COLUMN the group after it: the top row is what
    happens after a scatter out of group 0, the next row out of group 1, and so
    on. Every row is normalized on its own, so rows only need non-negative
    entries and a positive sum (different rows may have different totals).
    Nothing forces a direction: which entries are down- or up-scatter is
    decided entirely by the matrix. The row sums do NOT set how often a
    scatter happens -- that is still sigma_s[g]; the row only says where the
    neutron goes. A row of all zeros is accepted only for a group with
    sigma_s[g] == 0 (that group never scatters, so the row is never used).
    Left out (None or []), a scatter out of group g lands uniformly in groups
    0..g.

    Absorption is implicit: sigma_a[g] = sigma_t[g] - sigma_s[g] - sigma_f[g],
    for every group it must be >= 0 (sigma_s + sigma_f can't exceed sigma_t).
    """

    def __init__(self, name, sigma_t, sigma_s, sigma_f=None, nu=None, chi_fission=None, chi_scatter=None):
        name = canonical_material_name(name)  # raises ValueError early if unknown
        sigma_t = list(sigma_t)
        sigma_s = list(sigma_s)
        n_groups = len(sigma_t)

        if n_groups == 0:
            raise ValueError(f"{name}: sigma_t can't be empty")
        if len(sigma_s) != n_groups:
            raise ValueError(f"{name}: sigma_s must have {n_groups} elements (same as sigma_t), got {len(sigma_s)}")

        sigma_f = list(sigma_f) if sigma_f is not None else [0.0] * n_groups
        nu = list(nu) if nu is not None else [0.0] * n_groups
        if len(sigma_f) != n_groups:
            raise ValueError(f"{name}: sigma_f must have {n_groups} elements, got {len(sigma_f)}")
        if len(nu) != n_groups:
            raise ValueError(f"{name}: nu must have {n_groups} elements, got {len(nu)}")

        for g in range(n_groups):
            if sigma_t[g] < 0 or sigma_s[g] < 0 or sigma_f[g] < 0:
                raise ValueError(f"{name}: sigma values can't be negative (group {g})")
            if sigma_s[g] + sigma_f[g] > sigma_t[g] + 1e-12:
                raise ValueError(
                    f"{name}: sigma_s + sigma_f exceeds sigma_t at group {g} "
                    f"({sigma_s[g]} + {sigma_f[g]} > {sigma_t[g]})"
                )
            if nu[g] < 0:
                raise ValueError(f"{name}: nu can't be negative (group {g})")

        # chi_fission: None / [] = not given (children keep the parent's group).
        # The raw weights are kept as the caller typed them (so save/load round-trips
        # exactly); normalization happens in the C++ engine.
        if chi_fission is not None:
            chi_fission = [float(c) for c in chi_fission]
        if not chi_fission:
            chi_fission = None
        else:
            if len(chi_fission) != n_groups:
                raise ValueError(f"{name}: chi_fission must have {n_groups} elements (one per group), got {len(chi_fission)}")
            for g, c in enumerate(chi_fission):
                if not math.isfinite(c) or c < 0:
                    raise ValueError(f"{name}: chi_fission must be finite and non-negative (group {g}: {c})")
            if sum(chi_fission) <= 0:
                raise ValueError(f"{name}: chi_fission must have a positive sum (can't normalize all zeros)")

        # chi_scatter: None / [] = not given (uniform over groups 0..g). Raw rows are kept as
        # typed (so save/load round-trips exactly); normalization happens in the C++ engine.
        if chi_scatter is not None:
            try:
                chi_scatter = [[float(v) for v in row] for row in chi_scatter]
            except TypeError:
                raise ValueError(
                    f"{name}: chi_scatter must be an n_groups x n_groups matrix (a list of lists), "
                    f"e.g. [[0.5, 0.5], [0.0, 1.0]] for 2 groups"
                ) from None
        if not chi_scatter:
            chi_scatter = None
        else:
            if len(chi_scatter) != n_groups:
                raise ValueError(
                    f"{name}: chi_scatter must have {n_groups} rows (one per group before the scatter), got {len(chi_scatter)}"
                )
            for i, row in enumerate(chi_scatter):
                if len(row) != n_groups:
                    raise ValueError(f"{name}: chi_scatter row {i} must have {n_groups} elements, got {len(row)}")
                for j, v in enumerate(row):
                    if not math.isfinite(v) or v < 0:
                        raise ValueError(f"{name}: chi_scatter must be finite and non-negative (row {i}, column {j}: {v})")
                if sum(row) <= 0 and sigma_s[i] > 0:
                    raise ValueError(
                        f"{name}: chi_scatter row {i} sums to zero, but sigma_s[{i}] = {sigma_s[i]} > 0 "
                        f"(a neutron scattering out of group {i} would have nowhere to go)"
                    )

        self.name = name
        self.sigma_t = sigma_t
        self.sigma_s = sigma_s
        self.sigma_f = sigma_f
        self.nu = nu
        self.chi_fission = chi_fission
        self.chi_scatter = chi_scatter
        self.n_groups = n_groups


class MaterialLibrary:
    """Container for several GroupMaterial definitions, saved/loaded as ONE XML file.

        lib = MaterialLibrary()
        lib.add_material(name="Fe", sigma_t=[...], sigma_s=[...])
        lib.add_material(name="C",  sigma_t=[...], sigma_s=[...])
        lib.export("materials.xml")          # single file, every material inside

        lib2 = MaterialLibrary.load("materials.xml")

    Every material added must have the SAME number of groups as the ones
    already in the library (the whole system needs one consistent group
    structure) -- checked immediately in add_material(), not later at run().
    """

    def __init__(self):
        self.materials = {}  # name -> GroupMaterial
        self.n_groups = None

    def add_material(self, name, sigma_t, sigma_s, sigma_f=None, nu=None, chi_fission=None, chi_scatter=None):
        gm = GroupMaterial(name, sigma_t, sigma_s, sigma_f, nu, chi_fission, chi_scatter)
        if self.n_groups is None:
            self.n_groups = gm.n_groups
        elif gm.n_groups != self.n_groups:
            raise ValueError(
                f"'{name}' has {gm.n_groups} groups, but this library already uses "
                f"{self.n_groups} groups -- every material must share the same group structure"
            )
        self.materials[gm.name] = gm  # keyed by the CANONICAL name, not necessarily what the caller typed
        return self  # chainable

    def __contains__(self, name):
        try:
            name = canonical_material_name(name)
        except ValueError:
            return False  # not a material MCMR recognizes at all -> definitely not in here
        return name in self.materials

    def __getitem__(self, name):
        return self.materials[canonical_material_name(name)]

    # ------------------------------------------------------------------ #
    # Save / load to XML -- one file holds every material
    # ------------------------------------------------------------------ #
    def export(self, filename):
        root = ET.Element("mcmr_group_materials")
        root.set("n_groups", str(self.n_groups or 0))

        for gm in self.materials.values():
            m_el = ET.SubElement(root, "material")
            m_el.set("name", gm.name)
            ET.SubElement(m_el, "sigma_t").text = ",".join(map(str, gm.sigma_t))
            ET.SubElement(m_el, "sigma_s").text = ",".join(map(str, gm.sigma_s))
            ET.SubElement(m_el, "sigma_f").text = ",".join(map(str, gm.sigma_f))
            ET.SubElement(m_el, "nu").text = ",".join(map(str, gm.nu))
            if gm.chi_fission is not None:  # only written when given, so "absent" survives a round trip
                ET.SubElement(m_el, "chi_fission").text = ",".join(map(str, gm.chi_fission))
            if gm.chi_scatter is not None:  # one <row> per group before the scatter (row 0 first)
                cs_el = ET.SubElement(m_el, "chi_scatter")
                for row in gm.chi_scatter:
                    ET.SubElement(cs_el, "row").text = ",".join(map(str, row))

        tree = ET.ElementTree(root)
        ET.indent(tree, space="  ")
        tree.write(filename, xml_declaration=True, encoding="UTF-8")
        return filename

    @classmethod
    def load(cls, filename):
        tree = ET.parse(filename)
        root = tree.getroot()

        lib = cls()
        for m_el in root.findall("material"):
            name = m_el.get("name")
            sigma_t = [float(v) for v in m_el.find("sigma_t").text.split(",") if v]
            sigma_s = [float(v) for v in m_el.find("sigma_s").text.split(",") if v]
            sf_el = m_el.find("sigma_f")
            nu_el = m_el.find("nu")
            sigma_f = [float(v) for v in sf_el.text.split(",") if v] if sf_el is not None and sf_el.text else None
            nu = [float(v) for v in nu_el.text.split(",") if v] if nu_el is not None and nu_el.text else None
            chi_el = m_el.find("chi_fission")  # optional: files saved before chi_fission existed have none
            chi_fission = [float(v) for v in chi_el.text.split(",") if v] if chi_el is not None and chi_el.text else None
            cs_el = m_el.find("chi_scatter")  # optional: files saved before chi_scatter existed have none
            chi_scatter = None
            if cs_el is not None:
                rows = [[float(v) for v in r.text.split(",") if v] for r in cs_el.findall("row") if r.text]
                chi_scatter = rows or None
            lib.add_material(name=name, sigma_t=sigma_t, sigma_s=sigma_s, sigma_f=sigma_f, nu=nu,
                             chi_fission=chi_fission, chi_scatter=chi_scatter)

        return lib
