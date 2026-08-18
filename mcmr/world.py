import xml.etree.ElementTree as ET

__all__ = ["World"]


class World:
    """Representasi tunggal sebuah dunia simulasi -- geometri, grid, material, sumber, boundary.

    Ini adalah SATU-SATUNYA jalur, tidak peduli dunia itu didefinisikan manual
    (array mentah) atau lewat Geometry builder:

        world = mcmr.World(x_world=50, y_world=50, x_grid=[...], y_grid=[...],
                            material_matrix=[...], sources=[...])
        world = geom.build()                    # dari Geometry, hasilnya juga World

    Dari titik ini, jalurnya identik:

        sim = world.run(N=5000)                 # jalankan simulasi
        world.export("world.xml")               # simpan definisi dunia ke XML
        world2 = mcmr.World.load("world.xml")    # baca lagi -> World baru, bisa .run() lagi
    """

    def __init__(self, x_world, y_world, x_grid, y_grid, material_matrix, sources,
                 bc_top="vacuum", bc_bot="vacuum", bc_left="vacuum", bc_right="vacuum"):
        ny = len(material_matrix)
        nx = len(material_matrix[0]) if ny else 0
        if any(len(row) != nx for row in material_matrix):
            raise ValueError("setiap baris material_matrix harus punya jumlah kolom yang sama")
        if len(sources) != ny or any(len(row) != nx for row in sources):
            raise ValueError("sources harus punya shape yang sama persis dengan material_matrix")

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
    # Jalankan simulasi
    # ------------------------------------------------------------------ #
    def run(self, N, max_save=50):
        """Jalankan simulasi Monte Carlo untuk dunia ini. Return objek Simulation (sudah di-.run()).

        N        : jumlah partikel neutron yang disimulasikan
        max_save : jumlah maksimum lintasan neutron yang disimpan untuk plotting

        Konvensi index [row][col] material_matrix / sources: row=0 adalah baris
        PALING ATAS (y tertinggi), col=0 adalah kolom PALING KIRI (x=0) -- ditulis
        natural seperti menggambar grid di kertas.

        bc_top, bc_bot, bc_left, bc_right : "vacuum" (neutron mati/leak) atau
        "reflective" (neutron memantul balik, energi tetap).
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
    # Simpan / baca ke XML
    # ------------------------------------------------------------------ #
    def export(self, filename):
        """Simpan dunia ini ke file XML (mcmr_world). Return filename yang ditulis.

        File ini TERPISAH dari XML hasil simulasi (mcmr_results.xml) -- 1 file = 1 informasi.
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
        """Baca file XML hasil World.export(), kembalikan World baru."""
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