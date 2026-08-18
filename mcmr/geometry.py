"""Geometry builder -- a concise & flexible way to define a simulation world.

Instead of manually writing x_grid, y_grid, material_matrix, and sources
cell by cell, the user simply does:

    geom = Geometry(x_world=50, y_world=50, default_material="Fe", default_source=0)
    geom.add_region(0, 0, 15, 20, material="Pb", source=1)
    geom.fill_row(12, 15, material="Pb")          # x from 12..15, full y height
    geom.fill_col(10, 20, source=2)                # y from 10..20, full x width

The same command (add_region / fill_row / fill_col) is used to change
material only, source only, or both at once -- depending on which
parameter is filled in. A parameter left unset (None) is NOT changed from
its default / previous value. A region registered later overrides a region
that overlaps with an earlier one (painter's algorithm), just like
CSS/drawing layer over layer.

The grid (x_grid, y_grid) doesn't need to be defined manually -- it's
automatically inferred from all region boundary coordinates ever registered.

geom.build() produces a World object -- identical to a World built
manually. From there the path is identical: world.run(N), world.export(filename),
etc.; there's no separate method for "result from Geometry" vs "manual".
"""

from .world import World

__all__ = ["Geometry"]


class Geometry:
    def __init__(self, x_world, y_world, default_material,
                 default_source=0,
                 bc_top="vacuum", bc_bot="vacuum", bc_left="vacuum", bc_right="vacuum"):
        if x_world <= 0 or y_world <= 0:
            raise ValueError("x_world and y_world must be positive")

        self.x_world = float(x_world)
        self.y_world = float(y_world)
        self.default_material = default_material
        self.default_source = default_source
        self.bc_top = bc_top
        self.bc_bot = bc_bot
        self.bc_left = bc_left
        self.bc_right = bc_right

        self._ops = []  # list of dict(x1,y1,x2,y2,material,source), in call order
        self._x_bounds = {0.0, self.x_world}
        self._y_bounds = {0.0, self.y_world}

    # ------------------------------------------------------------------ #
    # Region commands
    # ------------------------------------------------------------------ #
    def add_region(self, x1, y1, x2, y2, material=None, source=None):
        """Define a rectangle from (x1, y1) to (x2, y2).

        material : material name (str) for this region. None = unchanged.
        source   : neutron source weight (number) for this region. None = unchanged.
        At least one of material/source must be provided.
        """
        x1, x2 = float(x1), float(x2)
        y1, y2 = float(y1), float(y2)

        if material is None and source is None:
            raise ValueError("add_region needs at least one of material or source")
        if not (0 <= x1 < x2 <= self.x_world):
            raise ValueError(f"x1={x1}, x2={x2} out of world bounds [0, {self.x_world}] or x1 >= x2")
        if not (0 <= y1 < y2 <= self.y_world):
            raise ValueError(f"y1={y1}, y2={y2} out of world bounds [0, {self.y_world}] or y1 >= y2")

        self._ops.append({"x1": x1, "y1": y1, "x2": x2, "y2": y2,
                           "material": material, "source": source})
        self._x_bounds.update((x1, x2))
        self._y_bounds.update((y1, y2))
        return self  # chainable: geom.add_region(...).add_region(...)

    def fill_row(self, x1, x2, material=None, source=None):
        """Fill a vertical strip: x from x1 to x2, covering the FULL y height."""
        return self.add_region(x1, 0, x2, self.y_world, material=material, source=source)

    def fill_col(self, y1, y2, material=None, source=None):
        """Fill a horizontal strip: y from y1 to y2, covering the FULL x width."""
        return self.add_region(0, y1, self.x_world, y2, material=material, source=source)

    # ------------------------------------------------------------------ #
    # Build -> World
    # ------------------------------------------------------------------ #
    def build(self):
        """Compile all registered regions into a World object.

        Returns a World -- the subsequent path (world.run(N), world.export(filename))
        is identical to a World built manually.
        """
        x_edges = sorted(self._x_bounds)
        y_edges = sorted(self._y_bounds)
        nx = len(x_edges) - 1
        ny = len(y_edges) - 1

        material_matrix = [[self.default_material] * nx for _ in range(ny)]
        sources = [[self.default_source] * nx for _ in range(ny)]

        for op in self._ops:
            for row in range(ny):
                # row=0 = topmost = highest y
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

        x_grid = x_edges[1:-1]  # drop 0 and x_world -- those are internal boundaries only
        y_grid = y_edges[1:-1]

        return World(
            x_world=self.x_world, y_world=self.y_world,
            x_grid=x_grid, y_grid=y_grid,
            material_matrix=material_matrix, sources=sources,
            bc_top=self.bc_top, bc_bot=self.bc_bot,
            bc_left=self.bc_left, bc_right=self.bc_right,
        )
