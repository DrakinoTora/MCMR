import xml.etree.ElementTree as ET
import matplotlib.pyplot as plt
import numpy as np

class ResultsPlotter:
    """Kelas untuk membaca XML luaran MCMR dan memvisualisasikan grafik."""
    def __init__(self, xml_filename="mcmr_results.xml"):
        self.xml_filename = xml_filename
        self.tree = ET.parse(xml_filename)
        self.root = self.tree.getroot()

    def plot_trajectories(self, box_bounds=None, world_bounds=None, cmap_name='plasma'):
        born_text = self.root.find(".//energy_born").text
        e_born_all = np.array(list(map(float, born_text.split(',')))) / 1e6  # eV -> MeV

        histories = self.root.findall(".//particle_history")

        e_born_plotted = np.array([e_born_all[int(h.get("id"))] for h in histories])

        cmap = plt.get_cmap(cmap_name)
        norm = plt.Normalize(vmin=e_born_plotted.min(), vmax=e_born_plotted.max())

        fig, ax = plt.subplots(figsize=(7, 6))
        for history, e in zip(histories, e_born_plotted):
            x_vals = list(map(float, history.find("x").text.split(',')))
            y_vals = list(map(float, history.find("y").text.split(',')))
            ax.plot(x_vals, y_vals, '-o', markersize=2, alpha=0.7, color=cmap(norm(e)))

        if box_bounds:
            x1, y1, x2, y2 = box_bounds
            ax.add_patch(plt.Rectangle((x1, y1), x2 - x1, y2 - y1,
                                        fill=None, edgecolor='red', linewidth=2, label='Box Material'))
        if world_bounds:
            ax.set_xlim(0, world_bounds[0])
            ax.set_ylim(0, world_bounds[1])

        sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
        sm.set_array([])
        fig.colorbar(sm, ax=ax, label='Energy Borm (MeV)')

        ax.set_xlabel('X (cm)')
        ax.set_ylabel('Y (cm)')
        ax.set_title('Neutron Path')
        ax.legend()
        ax.grid(True)
        plt.show()

    def plot_energy_born(self):
        text = self.root.find(".//energy_born").text
        if text:
            e_born = np.array(list(map(float, text.split(',')))) / 1e6 # MeV
            plt.figure(figsize=(7, 4))
            plt.hist(e_born, bins=50, color='orange', alpha=0.7, density=True)
            plt.xlabel('Energy (MeV)')
            plt.ylabel('Density Probability')
            plt.title('Born Energy')
            plt.grid(True)
            plt.show()

    def plot_energy_leak(self):
        text = self.root.find(".//energy_leak").text
        if text:
            e_leak = np.array(list(map(float, text.split(',')))) / 1e6 # MeV
            plt.figure(figsize=(7, 4))
            plt.hist(e_leak, bins=50, color='green', alpha=0.7)
            plt.xlabel('Energy (MeV)')
            plt.ylabel('Particle')
            plt.title('Leak Energy')
            plt.grid(True)
            plt.show()