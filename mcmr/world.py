import xml.etree.ElementTree as ET

__all__ = ["World"]


class World:
    """Single representation of a simulation world -- geometry, grid, materials, sources, boundary.

    This is the ONLY path, regardless of whether the world was defined manually
    (raw arrays) or via the Geometry builder:

        world = mcmr.World(x_world=50, y_world=50, x_grid=[...], y_grid=[...],
                            material_matrix=[...], sources=[...])
        world = geom.build()                    # from Geometry, the result is also a World

    From this point on, the path is identical:

        sim = world.run(N=5000)                 # run the simulation
        world.export("world.xml")               # save the world definition to XML
        world2 = mcmr.World.load("world.xml")    # read it back -> new World, can .run() again

    A world can also hold true circular regions, added on top of the rectangular
    grid via add_circle(). These are real circles in the physics engine (checked
    with an actual line-circle intersection, not a rasterized/pixelated approximation):

        world.add_circle(x=25, y=25, r=5, material="C", source=1)
    """

    def __init__(self, x_world, y_world, x_grid, y_grid, material_matrix, sources,
                 bc_top="vacuum", bc_bot="vacuum", bc_left="vacuum", bc_right="vacuum"):
        ny = len(material_matrix)
        nx = len(material_matrix[0]) if ny else 0
        if any(len(row) != nx for row in material_matrix):
            raise ValueError("every row of material_matrix must have the same number of columns")
        if len(sources) != ny or any(len(row) != nx for row in sources):
            raise ValueError("sources must have the exact same shape as material_matrix")

        self.x_world = x_world
        self.y_world = y_world
        self.x_grid = list(x_grid)
        self.y_grid = list(y_grid)
        self.material_matrix = material_matrix
        self.sources = sources
        self.bc_top = bc_top
        self.bc_bot = bc_bot
        self.bc_left = bc_left
        self.bc_right = bc_right
        self.circles = []  # list of dict(cx, cy, r, material, source)

    # ------------------------------------------------------------------ #
    # Circles -- true geometry, drawn on top of the rectangular grid
    # ------------------------------------------------------------------ #
    def add_circle(self, x, y, r, material, source=1):
        """Add a true circular region on top of this world.

        This is real circular geometry in the physics engine (an actual
        line-circle intersection is solved during transport) -- NOT a
        rasterized/pixelated approximation made of small rectangles.

        x, y     : center of the circle.
        r        : radius. Must fit entirely inside the world bounds, or a
                   ValueError is raised.
        material : material name for inside the circle.
        source   : neutron source weight for this circle (default 1). Use 0
                   if this circle should not emit any neutrons.

        Circles added later are drawn ON TOP of earlier circles wherever they
        overlap (painter's algorithm) -- same convention as Geometry regions.
        """
        if r <= 0:
            raise ValueError("circle radius must be positive")
        if not (0 <= x - r and x + r <= self.x_world and 0 <= y - r and y + r <= self.y_world):
            raise ValueError(
                f"circle at ({x}, {y}) with radius {r} goes outside world bounds "
                f"[0, {self.x_world}] x [0, {self.y_world}]"
            )
        if source < 0:
            raise ValueError("circle source can't be negative")

        self.circles.append({"cx": x, "cy": y, "r": r, "material": material, "source": source})
        return self  # chainable

    def add_circles(self, radius_matrix, material_matrix, source_matrix=None):
        """Add many circles at once, one per grid cell, using the SAME [row][col]
        convention as this world's own material_matrix/sources (row=0 = topmost,
        col=0 = leftmost). Each circle is automatically centered in the middle
        of its grid cell -- internally this just calls add_circle() per cell,
        same validation applies.

        radius_matrix   : must have the exact same shape as this world's
                           material_matrix. A cell value of 0 (or any falsy
                           value) means "no circle in this cell".
        material_matrix : material name for each active cell (ignored where
                           radius_matrix is 0).
        source_matrix   : neutron source weight for each active cell.
                           Optional -- defaults to 1 for every active cell.

        A circle's radius can't exceed half of its own cell's shortest side
        (otherwise it would stick out of the cell it's centered in) -- raises
        ValueError naming the offending [row][col] if it does.
        """
        ny = len(self.material_matrix)
        nx = len(self.material_matrix[0]) if ny else 0

        if len(radius_matrix) != ny or any(len(row) != nx for row in radius_matrix):
            raise ValueError("radius_matrix must have the exact same shape as this world's material_matrix")
        if len(material_matrix) != ny or any(len(row) != nx for row in material_matrix):
            raise ValueError("material_matrix must have the exact same shape as this world's material_matrix")
        if source_matrix is not None and (
            len(source_matrix) != ny or any(len(row) != nx for row in source_matrix)
        ):
            raise ValueError("source_matrix must have the exact same shape as this world's material_matrix")

        x_edges = [0.0] + list(self.x_grid) + [self.x_world]
        y_edges = [0.0] + list(self.y_grid) + [self.y_world]

        for row in range(ny):
            cy1, cy2 = y_edges[ny - 1 - row], y_edges[ny - row]  # row=0 -> topmost -> highest y
            cy = (cy1 + cy2) / 2
            cell_h = cy2 - cy1
            for col in range(nx):
                r = radius_matrix[row][col]
                if not r:
                    continue

                cx1, cx2 = x_edges[col], x_edges[col + 1]
                cx = (cx1 + cx2) / 2
                cell_w = cx2 - cx1

                max_r = min(cell_w, cell_h) / 2
                if r > max_r:
                    raise ValueError(
                        f"radius {r} at cell [row={row}][col={col}] exceeds half of its cell's "
                        f"shortest side ({max_r}) -- the circle would stick out of its own cell"
                    )

                source = source_matrix[row][col] if source_matrix is not None else 1
                self.add_circle(x=cx, y=cy, r=r, material=material_matrix[row][col], source=source)

        return self

    # ------------------------------------------------------------------ #
    # Run simulation
    # ------------------------------------------------------------------ #
    def run(self, N, max_save=50, mode="continuous", group_materials=None):
        """Run the Monte Carlo simulation for this world. Returns a Simulation
        (or SimulationMG) object, already .run().

        N        : number of neutron particles to simulate
        max_save : maximum number of neutron trajectories saved for plotting
        mode     : "continuous" (default, unchanged behavior -- Watt spectrum
                   birth + interpolated continuous cross sections) or "group"
                   (discretized energy -- see group_materials below).
        group_materials : required when mode="group". A mcmr.MaterialLibrary
                   holding a mcmr.GroupMaterial for every material name used
                   in this world's material_matrix / circles.

        [row][col] index convention for material_matrix / sources: row=0 is the
        TOPMOST row (highest y), col=0 is the LEFTMOST column (x=0) -- written
        naturally, like drawing a grid on paper.

        bc_top, bc_bot, bc_left, bc_right : "vacuum" (neutron dies/leaks) or
        "reflective" (neutron bounces back, energy unchanged).
        """
        if mode not in ("continuous", "group"):
            raise ValueError(f"mode must be 'continuous' or 'group', got {mode!r}")

        if mode == "continuous":
            return self._run_continuous(N, max_save)
        return self._run_group(N, max_save, group_materials)

    def _run_continuous(self, N, max_save):
        from ._mcmr_cpp import Simulation
        from .cross_section import load_all_materials

        sim = Simulation(
            N=N, x_world=self.x_world, y_world=self.y_world,
            x_grid=self.x_grid, y_grid=self.y_grid,
            material_matrix=self.material_matrix, sources=self.sources,
            circle_cx=[c["cx"] for c in self.circles],
            circle_cy=[c["cy"] for c in self.circles],
            circle_r=[c["r"] for c in self.circles],
            circle_material=[c["material"] for c in self.circles],
            circle_source=[c["source"] for c in self.circles],
            max_history_save=max_save,
            bc_top=self.bc_top, bc_bot=self.bc_bot, bc_left=self.bc_left, bc_right=self.bc_right,
        )
        E_tot, Sig_tot, E_scat, Sig_scat = load_all_materials()
        sim.set_cross_sections(E_tot, Sig_tot, E_scat, Sig_scat)

        sim.run()
        return sim

    def _run_group(self, N, max_save, group_materials):
        from ._mcmr_cpp import SimulationMG, canonical_material_name

        if group_materials is None:
            raise ValueError("mode='group' requires group_materials=<a MaterialLibrary>")

        # Resolve every raw name actually used in this world (material_matrix + circles)
        # to the SAME canonical spelling the C++ engine keys its group data on
        # (via canonical_material_name(), a thin wrapper around get_material_info()).
        # This is what keeps material_matrix=[["fe"]] and
        # group_materials.add_material("Fe", ...) lined up even though they were
        # spelled differently -- and turns an unrecognized name into a clear
        # Python ValueError right here, instead of a std::map::at crash deep
        # inside the C++ transport loop once the simulation is already running.
        raw_names = {name for row in self.material_matrix for name in row}
        raw_names |= {c["material"] for c in self.circles}
        canonical_of = {}
        unknown = []
        for raw in raw_names:
            try:
                canonical_of[raw] = canonical_material_name(raw)
            except ValueError:
                unknown.append(raw)
        if unknown:
            raise ValueError(
                f"material_matrix/circles use unrecognized material name(s): {', '.join(sorted(unknown))} "
                "(known: Be, C, Fe, Pb, and their aliases -- see get_material_info() in src/material.cpp)"
            )

        used_names = set(canonical_of.values())
        missing = [n for n in used_names if n not in group_materials]
        if missing:
            raise ValueError(f"group_materials is missing data for: {', '.join(sorted(missing))}")

        sigma_t, sigma_s, sigma_f, nu = {}, {}, {}, {}
        for name in used_names:
            gm = group_materials[name]
            sigma_t[name] = gm.sigma_t
            sigma_s[name] = gm.sigma_s
            sigma_f[name] = gm.sigma_f
            nu[name] = gm.nu

        # material_matrix/circle_material go to C++ as-is (any alias); the C++
        # constructor resolves them to the same canonical symbols via
        # get_material_info(), so they land on the sigma_t/sigma_s/... keys above.
        sim = SimulationMG(
            N=N, x_world=self.x_world, y_world=self.y_world,
            x_grid=self.x_grid, y_grid=self.y_grid,
            material_matrix=self.material_matrix, sources=self.sources,
            circle_cx=[c["cx"] for c in self.circles],
            circle_cy=[c["cy"] for c in self.circles],
            circle_r=[c["r"] for c in self.circles],
            circle_material=[c["material"] for c in self.circles],
            circle_source=[c["source"] for c in self.circles],
            max_history_save=max_save,
            bc_top=self.bc_top, bc_bot=self.bc_bot, bc_left=self.bc_left, bc_right=self.bc_right,
        )
        sim.set_group_data(group_materials.n_groups, sigma_t, sigma_s, sigma_f, nu)
        sim.run()
        return sim

    # ------------------------------------------------------------------ #
    # Save / load to XML
    # ------------------------------------------------------------------ #
    def export(self, filename):
        """Save this world to an XML file (mcmr_world). Returns the filename written.

        This file is SEPARATE from the simulation results XML (mcmr_results.xml) -- 1 file = 1 piece of information.
        """
        root = ET.Element("mcmr_world")

        dims = ET.SubElement(root, "dimensions")
        dims.set("x_world", str(self.x_world))
        dims.set("y_world", str(self.y_world))

        grid_el = ET.SubElement(root, "grid")
        ET.SubElement(grid_el, "x_grid").text = ",".join(map(str, self.x_grid))
        ET.SubElement(grid_el, "y_grid").text = ",".join(map(str, self.y_grid))

        bc_el = ET.SubElement(root, "boundary")
        bc_el.set("top", self.bc_top)
        bc_el.set("bottom", self.bc_bot)
        bc_el.set("left", self.bc_left)
        bc_el.set("right", self.bc_right)

        mats_el = ET.SubElement(root, "materials")
        for i, row in enumerate(self.material_matrix):
            row_el = ET.SubElement(mats_el, "row")
            row_el.set("index", str(i))
            row_el.text = ",".join(row)

        src_el = ET.SubElement(root, "sources")
        for i, row in enumerate(self.sources):
            row_el = ET.SubElement(src_el, "row")
            row_el.set("index", str(i))
            row_el.text = ",".join(map(str, row))

        circles_el = ET.SubElement(root, "circles")
        for c in self.circles:
            c_el = ET.SubElement(circles_el, "circle")
            c_el.set("cx", str(c["cx"]))
            c_el.set("cy", str(c["cy"]))
            c_el.set("r", str(c["r"]))
            c_el.set("material", c["material"])
            c_el.set("source", str(c["source"]))

        tree = ET.ElementTree(root)
        ET.indent(tree, space="  ")
        tree.write(filename, xml_declaration=True, encoding="UTF-8")
        return filename

    @classmethod
    def load(cls, filename):
        """Read the XML file produced by World.export(), return a new World."""
        tree = ET.parse(filename)
        root = tree.getroot()

        dims = root.find("dimensions")
        x_world = float(dims.get("x_world"))
        y_world = float(dims.get("y_world"))

        grid_el = root.find("grid")
        x_grid = [float(v) for v in (grid_el.find("x_grid").text or "").split(",") if v]
        y_grid = [float(v) for v in (grid_el.find("y_grid").text or "").split(",") if v]

        bc_el = root.find("boundary")
        bc_top = bc_el.get("top", "vacuum")
        bc_bot = bc_el.get("bottom", "vacuum")
        bc_left = bc_el.get("left", "vacuum")
        bc_right = bc_el.get("right", "vacuum")

        material_matrix = [row.text.split(",") for row in root.find("materials").findall("row")]
        sources = [[float(v) for v in row.text.split(",")] for row in root.find("sources").findall("row")]

        world = cls(
            x_world=x_world, y_world=y_world,
            x_grid=x_grid, y_grid=y_grid,
            material_matrix=material_matrix, sources=sources,
            bc_top=bc_top, bc_bot=bc_bot, bc_left=bc_left, bc_right=bc_right,
        )

        circles_el = root.find("circles")
        if circles_el is not None:
            for c_el in circles_el.findall("circle"):
                world.add_circle(
                    x=float(c_el.get("cx")), y=float(c_el.get("cy")), r=float(c_el.get("r")),
                    material=c_el.get("material"), source=float(c_el.get("source")),
                )

        return world