import xml.etree.ElementTree as ET
import matplotlib.pyplot as plt
import numpy as np

class ResultsPlotter:
    """Kelas untuk membaca XML luaran MCMR dan memvisualisasikan grafik."""
    def __init__(self, xml_filename="mcmr_results.xml"):
        self.xml_filename = xml_filename
        self.tree = ET.parse(xml_filename)
        self.root = self.tree.getroot()

    def plot_trajectories(self, box_bounds=None, world_bounds=None):
        plt.figure(figsize=(7, 6))
        for history in self.root.findall(".//particle_history"):
            x_vals = list(map(float, history.find("x").text.split(',')))
            y_vals = list(map(float, history.find("y").text.split(',')))
            plt.plot(x_vals, y_vals, '-o', markersize=2, alpha=0.5)

        if box_bounds:
            x1, y1, x2, y2 = box_bounds
            plt.gca().add_patch(plt.Rectangle((x1, y1), x2 - x1, y2 - y1, 
                                              fill=None, edgecolor='red', linewidth=2, label='Box Material'))
        if world_bounds:
            plt.xlim(0, world_bounds[0])
            plt.ylim(0, world_bounds[1])

        plt.xlabel('X (cm)')
        plt.ylabel('Y (cm)')
        plt.title('Lintasan Neutron (Dari XML Hasil MCMR)')
        plt.legend()
        plt.grid(True)
        plt.show()

    def plot_energy_born(self):
        text = self.root.find(".//energy_born").text
        if text:
            e_born = np.array(list(map(float, text.split(',')))) / 1e6 # MeV
            plt.figure(figsize=(7, 4))
            plt.hist(e_born, bins=50, color='orange', alpha=0.7, density=True)
            plt.xlabel('Energi (MeV)')
            plt.ylabel('Densitas Probabilitas')
            plt.title('Distribusi Energi Neutron Dilahirkan')
            plt.grid(True)
            plt.show()

    def plot_energy_leak(self):
        text = self.root.find(".//energy_leak").text
        if text:
            e_leak = np.array(list(map(float, text.split(',')))) / 1e6 # MeV
            plt.figure(figsize=(7, 4))
            plt.hist(e_leak, bins=50, color='green', alpha=0.7)
            plt.xlabel('Energi (MeV)')
            plt.ylabel('Jumlah Partikel')
            plt.title('Distribusi Energi Neutron Bocor')
            plt.grid(True)
            plt.show()