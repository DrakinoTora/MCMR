import itertools
import xml.etree.ElementTree as ET
import numpy as np

import matplotlib.patheffects as pe
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Circle

from .world import World

MATERIAL_COLORS = {
    "Fe": "#cfd8dc",
    "Pb": "#bcaaa4",
    "C":  "#90a4ae",
    "Be": "#c8e6c9", 
}
_FALLBACK_PALETTE = itertools.cycle(plt.get_cmap("Pastel1").colors)

TRAJECTORY_PALETTE = [
    "#e6194b", "#3cb44b", "#4363d8", "#f58231", "#911eb4",
    "#f032e6", "#008080", "#9a6324", "#800000", "#000075",
]


class ResultsPlotter:
    """Visualization of neutron trajectories from MCMR simulation results."""

    def __init__(self, xml_filename="mcmr_results.xml"):
        self.xml_filename = xml_filename
        self.tree = ET.parse(xml_filename)
        self.root = self.tree.getroot()

    def _material_color(self, name):
        if name not in MATERIAL_COLORS:
            MATERIAL_COLORS[name] = next(_FALLBACK_PALETTE)
        return MATERIAL_COLORS[name]

    def plot_trajectories(self, world_xml=None, x_world=None, y_world=None,
                           x_grid=None, y_grid=None, material_matrix=None, circles=None):
        """
        Option 1 (recommended) -- read the grid from the world XML file produced by
        World.export(), no need to retype any parameters:
            plotter.plot_trajectories(world_xml="world.xml")

        Option 2 (manual) -- pass the grid parameters directly. The grid parameters
        must be EXACTLY the same as those used in world.run():
            plotter.plot_trajectories(x_world=..., y_world=..., x_grid=..., y_grid=..., material_matrix=..., circles=...)

        material_matrix : material_matrix[row][col], row=0 is the TOPMOST row (highest y),
                           col=0 is the LEFTMOST column (x=0) -- same as writing a grid on paper.
        circles         : list of dict(cx, cy, r, material, ...), same shape as World.circles.
                           Drawn on top of the rectangular grid, same as the physics engine.
                           Only used in Option 2 -- Option 1 reads it from world_xml automatically.
        """
        if world_xml is not None:
            w = World.load(world_xml)
            x_world, y_world = w.x_world, w.y_world
            x_grid, y_grid = w.x_grid, w.y_grid
            material_matrix = w.material_matrix
            circles = w.circles
        elif None in (x_world, y_world, x_grid, y_grid, material_matrix):
            raise ValueError(
                "give world_xml=... (recommended) OR all of "
                "x_world, y_world, x_grid, y_grid, material_matrix manually"
            )
        circles = circles or []

        x_edges = [0.0] + list(x_grid) + [x_world]
        y_edges = [0.0] + list(y_grid) + [y_world]
        nx = len(x_edges) - 1
        ny = len(y_edges) - 1

        fig, ax = plt.subplots(figsize=(8, 7))

        # background (rectangular grid)
        drawn_materials = {}
        for row in range(ny):
            iy = ny - 1 - row
            for col in range(nx):
                ix = col
                mat = material_matrix[row][col]
                color = self._material_color(mat)
                ax.add_patch(Rectangle(
                    (x_edges[ix], y_edges[iy]),
                    x_edges[ix + 1] - x_edges[ix],
                    y_edges[iy + 1] - y_edges[iy],
                    facecolor=color, zorder=0,
                ))
                drawn_materials[mat] = color

        # circular regions -- drawn on top of the rectangular grid, same as the
        # physics engine (painter's algorithm: later circles override earlier ones)
        for c in circles:
            mat = c["material"]
            color = self._material_color(mat)
            ax.add_patch(Circle(
                (c["cx"], c["cy"]), c["r"],
                facecolor=color, edgecolor="black", linewidth=0.6, zorder=1,
            ))
            drawn_materials[mat] = color

        # neutron track
        histories = self.root.findall(".//particle_history")
        for idx, history in enumerate(histories):
            x_vals = list(map(float, history.find("x").text.split(",")))
            y_vals = list(map(float, history.find("y").text.split(",")))
            color = TRAJECTORY_PALETTE[idx % len(TRAJECTORY_PALETTE)]
            ax.plot(
                x_vals, y_vals, "-o", markersize=3, linewidth=1.6,
                color=color, alpha=0.95, zorder=2,
                path_effects=[pe.withStroke(linewidth=3, foreground="white")],
            )

        # legend
        legend_handles = [
            Rectangle((0, 0), 1, 1, facecolor=c, edgecolor="gray", label=m)
            for m, c in drawn_materials.items()
        ]
        ax.legend(handles=legend_handles, title="Material",
                  loc="upper left", bbox_to_anchor=(1.02, 1.0))

        ax.set_xlim(0, x_world)
        ax.set_ylim(0, y_world)
        ax.set_xlabel("X (cm)")
        ax.set_ylabel("Y (cm)")
        ax.set_title("neutron trajectory")
        ax.set_aspect("equal")
        plt.tight_layout()
        plt.show()
    
    def _region_heatmap(self, tag, title, cmap, world_xml=None):
        el = self.root.find(f".//{tag}")
        if el is None:
            raise ValueError(
                f"this results XML has no <{tag}> block -- it was produced by an "
                "older engine build that didn't separate fission/absorption tallies; "
                "rerun the simulation with the current mcmr build."
            )

        nx = int(el.get("nx"))
        ny = int(el.get("ny"))
        dx = float(el.get("dx"))
        dy = float(el.get("dy"))
        rows = el.findall("row")
        # <row index="0"> = BOTTOM row (y=0), same as imshow(origin="lower")
        data = np.array([[int(v) for v in row.text.split(",")] for row in rows])

        x_world = nx * dx
        y_world = ny * dy

        fig, ax = plt.subplots(figsize=(8, 7))
        im = ax.imshow(
            data, origin="lower", extent=[0, x_world, 0, y_world],
            cmap=cmap, interpolation="nearest", zorder=0,
        )
        cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
        cbar.set_label("Event count")

        if world_xml is not None:
            w = World.load(world_xml)
            x_edges = [0.0] + list(w.x_grid) + [w.x_world]
            y_edges = [0.0] + list(w.y_grid) + [w.y_world]
            for xe in x_edges:
                ax.axvline(xe, color="white", linewidth=0.5, alpha=0.4, zorder=1)
            for ye in y_edges:
                ax.axhline(ye, color="white", linewidth=0.5, alpha=0.4, zorder=1)
            for c in w.circles:
                ax.add_patch(Circle(
                    (c["cx"], c["cy"]), c["r"],
                    fill=False, edgecolor="white", linewidth=0.8, alpha=0.6, zorder=1,
                ))
            x_world, y_world = w.x_world, w.y_world

        if data.max() == 0:
            print(f"[mcmr] Note: <{tag}> is entirely zero in this results XML -- "
                "either no such event occurred, or (for fission_tally on a "
                "continuous-mode run) the engine has no fission model at all "
                "in that mode, so it's expected to stay zero.")

        ax.set_xlim(0, x_world)
        ax.set_ylim(0, y_world)
        ax.set_xlabel("X (cm)")
        ax.set_ylabel("Y (cm)")
        ax.set_title(title)
        ax.set_aspect("equal")
        plt.tight_layout()
        plt.show()

    def plot_fission_heatmap(self, world_xml=None):
        """Heatmap of where fission events happened (<fission_tally>).
        Group mode only; always all-zero for continuous mode.
        world_xml: optional, overlays grid lines + circle borders."""
        self._region_heatmap("fission_tally", "fission event heatmap", "hot", world_xml=world_xml)

    def plot_absorption_heatmap(self, world_xml=None):
        """Heatmap of where (non-fission) absorption happened (<absorp_tally>).
        Works for both continuous and group mode."""
        self._region_heatmap("absorp_tally", "absorption event heatmap", "viridis", world_xml=world_xml)

    def plot_k_generations(self, skip=None):
        """Plot k per generation (<k_generations>, group mode only).

        Draws k_estimate vs generation number, with a k = 1 (critical) reference
        line and the initial k the first generation divided nu by (k_used of
        generation 1) as a hollow marker at generation 0.

        skip : optional int. If given, the first `skip` generations are treated as
                inactive/warm-up and a dashed line shows the mean k_estimate of the
                remaining generations (with the mean +/- std in the legend).
        """
        el = self.root.find(".//k_generations")
        gens = el.findall("generation") if el is not None else []
        if not gens:
            print("[mcmr] Note: <k_generations> is empty or missing in this results XML "
                    "-- k is only recorded by group-mode runs, so there is nothing to plot.")
            return

        ids = [int(g.get("id")) for g in gens]
        k_used = [float(g.get("k_used")) for g in gens]
        k_est = [float(g.get("k_estimate")) for g in gens]

        fig, ax = plt.subplots(figsize=(8, 5))
        ax.axhline(1.0, color="gray", linestyle=":", linewidth=1.2)
        ax.plot(ids, k_est, "-o", markersize=4, linewidth=1.6, color="#4363d8")
        #ax.plot([0], [k_used[0]], "o", markerfacecolor="none", color="#4363d8")

        if skip is not None and 0 <= skip < len(k_est):
            active = np.array(k_est[skip:])
            mean, std = active.mean(), active.std()
            ax.axhline(mean, color="#e6194b", linestyle="--", linewidth=1.2)
            ax.axvspan(0, skip + 0.5, color="gray", alpha=0.12, zorder=0)

        ax.set_xlim(left=-0.5)
        ax.set_xlabel("Generation")
        ax.set_ylabel("k")
        ax.set_title(f"k per generation (final k = {k_est[-1]:.4f})")
        ax.grid(alpha=0.3)
        ax.legend()
        plt.tight_layout()
        plt.show()