"""Geometry builder -- cara ringkas & fleksibel untuk mendefinisikan dunia simulasi.

Daripada menulis manual x_grid, y_grid, material_matrix, dan sources sel-per-sel,
user cukup:

    geom = Geometry(x_world=50, y_world=50, default_material="Fe", default_source=0)
    geom.add_region(0, 0, 15, 20, material="Pb", source=1)
    geom.fill_row(12, 15, material="Pb")          # x dari 12..15, seluruh tinggi y
    geom.fill_col(10, 20, source=2)                # y dari 10..20, seluruh lebar x

Satu command yang sama (add_region / fill_row / fill_col) dipakai untuk mengubah
material saja, source saja, atau keduanya sekaligus -- tergantung parameter mana
yang diisi. Parameter yang tidak diisi (None) TIDAK diubah dari nilai default /
nilai sebelumnya. Region yang didaftarkan belakangan menimpa (override) region
yang tumpang tindih dengan region sebelumnya (painter's algorithm), sama seperti
CSS/menggambar layer di atas layer.

Grid (x_grid, y_grid) tidak perlu didefinisikan manual -- otomatis di-infer dari
semua koordinat batas region yang pernah didaftarkan.

geom.build() menghasilkan objek World -- sama persis dengan World yang dibuat
manual. Dari situ jalurnya identik: world.run(N), world.export(filename), dst,
tidak ada method terpisah untuk "hasil dari Geometry" vs "manual".
"""

from .world import World

__all__ = ["Geometry"]


class Geometry:
    def __init__(self, x_world, y_world, default_material,
                 default_source=0,
                 bc_top="vacuum", bc_bot="vacuum", bc_left="vacuum", bc_right="vacuum"):
        if x_world <= 0 or y_world <= 0:
            raise ValueError("x_world dan y_world harus positif")

        self.x_world = float(x_world)
        self.y_world = float(y_world)
        self.default_material = default_material
        self.default_source = default_source
        self.bc_top = bc_top
        self.bc_bot = bc_bot
        self.bc_left = bc_left
        self.bc_right = bc_right

        self._ops = []  # list of dict(x1,y1,x2,y2,material,source), urut sesuai pemanggilan
        self._x_bounds = {0.0, self.x_world}
        self._y_bounds = {0.0, self.y_world}

    # ------------------------------------------------------------------ #
    # Region commands
    # ------------------------------------------------------------------ #
    def add_region(self, x1, y1, x2, y2, material=None, source=None):
        """Definisikan rectangle dari (x1, y1) hingga (x2, y2).

        material : nama material (str) untuk region ini. None = tidak diubah.
        source   : bobot sumber neutron (angka) untuk region ini. None = tidak diubah.
        Minimal salah satu dari material/source harus diisi.
        """
        x1, x2 = float(x1), float(x2)
        y1, y2 = float(y1), float(y2)

        if material is None and source is None:
            raise ValueError("add_region butuh setidaknya salah satu dari material atau source")
        if not (0 <= x1 < x2 <= self.x_world):
            raise ValueError(f"x1={x1}, x2={x2} di luar batas dunia [0, {self.x_world}] atau x1 >= x2")
        if not (0 <= y1 < y2 <= self.y_world):
            raise ValueError(f"y1={y1}, y2={y2} di luar batas dunia [0, {self.y_world}] atau y1 >= y2")

        self._ops.append({"x1": x1, "y1": y1, "x2": x2, "y2": y2,
                           "material": material, "source": source})
        self._x_bounds.update((x1, x2))
        self._y_bounds.update((y1, y2))
        return self  # chainable: geom.add_region(...).add_region(...)

    def fill_row(self, x1, x2, material=None, source=None):
        """Isi strip vertikal: x dari x1 hingga x2, mencakup SELURUH tinggi y."""
        return self.add_region(x1, 0, x2, self.y_world, material=material, source=source)

    def fill_col(self, y1, y2, material=None, source=None):
        """Isi strip horizontal: y dari y1 hingga y2, mencakup SELURUH lebar x."""
        return self.add_region(0, y1, self.x_world, y2, material=material, source=source)

    # ------------------------------------------------------------------ #
    # Build -> World
    # ------------------------------------------------------------------ #
    def build(self):
        """Compile semua region yang terdaftar menjadi objek World.

        Return World -- jalur selanjutnya (world.run(N), world.export(filename))
        identik dengan World yang dibuat manual.
        """
        x_edges = sorted(self._x_bounds)
        y_edges = sorted(self._y_bounds)
        nx = len(x_edges) - 1
        ny = len(y_edges) - 1

        material_matrix = [[self.default_material] * nx for _ in range(ny)]
        sources = [[self.default_source] * nx for _ in range(ny)]

        for op in self._ops:
            for row in range(ny):
                # row=0 = paling atas = y tertinggi
                cy1, cy2 = y_edges[ny - 1 - row], y_edges[ny - row]
                cy_center = (cy1 + cy2) / 2
                if not (op["y1"] <= cy_center <= op["y2"]):
                    continue
                for col in range(nx):
                    cx1, cx2 = x_edges[col], x_edges[col + 1]
                    cx_center = (cx1 + cx2) / 2
                    if not (op["x1"] <= cx_center <= op["x2"]):
                        continue
                    if op["material"] is not None:
                        material_matrix[row][col] = op["material"]
                    if op["source"] is not None:
                        sources[row][col] = op["source"]

        x_grid = x_edges[1:-1]  # buang 0 dan x_world -- itu batas internal saja
        y_grid = y_edges[1:-1]

        return World(
            x_world=self.x_world, y_world=self.y_world,
            x_grid=x_grid, y_grid=y_grid,
            material_matrix=material_matrix, sources=sources,
            bc_top=self.bc_top, bc_bot=self.bc_bot,
            bc_left=self.bc_left, bc_right=self.bc_right,
        )