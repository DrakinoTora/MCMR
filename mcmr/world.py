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

    # ------------------------------------------------------------------ #
    # Run simulation
    # ------------------------------------------------------------------ #
    def run(self, N, max_save=50):
        """Run the Monte Carlo simulation for this world. Returns a Simulation object (already .run()).

        N        : number of neutron particles to simulate
        max_save : maximum number of neutron trajectories saved for plotting

        [row][col] index convention for material_matrix / sources: row=0 is the
        TOPMOST row (highest y), col=0 is the LEFTMOST column (x=0) -- written
        naturally, like drawing a grid on paper.

        bc_top, bc_bot, bc_left, bc_right : "vacuum" (neutron dies/leaks) or
        "reflective" (neutron bounces back, energy unchanged).
        """
        from ._mcmr_cpp import Simulation
        from .cross_section import load_all_materials

        sim = Simulation(N, self.x_world, self.y_world, self.x_grid, self.y_grid,
                          self.material_matrix, self.sources, max_save,
                          self.bc_top, self.bc_bot, self.bc_left, self.bc_right)
        E_tot, Sig_tot, E_scat, Sig_scat = load_all_materials()
        sim.set_cross_sections(E_tot, Sig_tot, E_scat, Sig_scat)

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
        x_grid = [float(v) for v in grid_el.find("x_grid").text.split(",") if v]
        y_grid = [float(v) for v in grid_el.find("y_grid").text.split(",") if v]

        bc_el = root.find("boundary")
        bc_top = bc_el.get("top", "vacuum")
        bc_bot = bc_el.get("bottom", "vacuum")
        bc_left = bc_el.get("left", "vacuum")
        bc_right = bc_el.get("right", "vacuum")

        material_matrix = [row.text.split(",") for row in root.find("materials").findall("row")]
        sources = [[float(v) for v in row.text.split(",")] for row in root.find("sources").findall("row")]

        return cls(
            x_world=x_world, y_world=y_world,
            x_grid=x_grid, y_grid=y_grid,
            material_matrix=material_matrix, sources=sources,
            bc_top=bc_top, bc_bot=bc_bot, bc_left=bc_left, bc_right=bc_right,
        )