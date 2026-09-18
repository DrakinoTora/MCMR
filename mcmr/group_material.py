import xml.etree.ElementTree as ET

__all__ = ["GroupMaterial", "MaterialLibrary"]


class GroupMaterial:
    """Multi-group (discretized energy) cross-section data for ONE material.

    `name` must match EXACTLY the material symbol used in this world's
    material_matrix / circles -- "Be", "C", "Fe" or "Pb" (same canonical
    spelling the C++ engine reports, see get_material_info() in
    src/material.cpp). No alias table on the Python side: the engine already
    resolves aliases (e.g. "besi", "iron") down to one of these 4 canonical
    symbols internally, so as long as material_matrix and MaterialLibrary
    agree on the canonical spelling, the two line up with no extra layer.

    Number of groups (n_groups) is inferred from len(sigma_t) -- there is no
    separate "how many groups" parameter to keep in sync by hand.

        Fe = GroupMaterial(
            name="Fe",
            sigma_t=[0.5, 0.2, 0.1],
            sigma_s=[0.3, 0.1, 0.06],
            sigma_f=[0.001, 0.0005, 0.00001],   # optional, default 0 (no fission)
            nu=[2.7, 2.5, 2.3],                 # optional, default 0
        )

    Absorption is implicit: sigma_a[g] = sigma_t[g] - sigma_s[g] - sigma_f[g],
    for every group it must be >= 0 (sigma_s + sigma_f can't exceed sigma_t).
    """

    def __init__(self, name, sigma_t, sigma_s, sigma_f=None, nu=None):
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

        self.name = name
        self.sigma_t = sigma_t
        self.sigma_s = sigma_s
        self.sigma_f = sigma_f
        self.nu = nu
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

    def add_material(self, name, sigma_t, sigma_s, sigma_f=None, nu=None):
        gm = GroupMaterial(name, sigma_t, sigma_s, sigma_f, nu)
        if self.n_groups is None:
            self.n_groups = gm.n_groups
        elif gm.n_groups != self.n_groups:
            raise ValueError(
                f"'{name}' has {gm.n_groups} groups, but this library already uses "
                f"{self.n_groups} groups -- every material must share the same group structure"
            )
        self.materials[name] = gm
        return self  # chainable

    def __contains__(self, name):
        return name in self.materials

    def __getitem__(self, name):
        return self.materials[name]

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
            lib.add_material(name=name, sigma_t=sigma_t, sigma_s=sigma_s, sigma_f=sigma_f, nu=nu)

        return lib
