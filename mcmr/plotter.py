import itertools
import xml.etree.ElementTree as ET

import matplotlib.patheffects as pe
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

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
    """Visualisasi lintasan neutron dari hasil simulasi MCMR."""

    def __init__(self, xml_filename="mcmr_results.xml"):
        self.xml_filename = xml_filename
        self.tree = ET.parse(xml_filename)
        self.root = self.tree.getroot()

    def _material_color(self, name):
        if name not in MATERIAL_COLORS:
            MATERIAL_COLORS[name] = next(_FALLBACK_PALETTE)
        return MATERIAL_COLORS[name]

    def plot_trajectories(self, x_world, y_world, x_grid, y_grid, material_matrix):
        """
        Gambar lintasan neutron di atas peta material grid (kotak-kotak warna warni).

        Parameter grid harus SAMA PERSIS dengan yang dipakai saat mcmr.load_and_run(),
        karena XML hasil simulasi tidak menyimpan info grid/material -- hanya lintasan & energi.

        x_world, y_world : ukuran dunia simulasi
        x_grid, y_grid   : batas interior grid (terurut naik)
        material_matrix  : material_matrix[i][j], i = indeks sel-x, j = indeks sel-y
        """
        x_edges = [0.0] + list(x_grid) + [x_world]
        y_edges = [0.0] + list(y_grid) + [y_world]

        fig, ax = plt.subplots(figsize=(8, 7))

        # background
        drawn_materials = {}
        for i in range(len(x_edges) - 1):
            for j in range(len(y_edges) - 1):
                mat = material_matrix[i][j]
                color = self._material_color(mat)
                ax.add_patch(Rectangle(
                    (x_edges[i], y_edges[j]),
                    x_edges[i + 1] - x_edges[i],
                    y_edges[j + 1] - y_edges[j],
                    facecolor=color, zorder=0,
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
