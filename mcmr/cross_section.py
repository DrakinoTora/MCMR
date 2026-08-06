import pandas as pd
from bs4 import BeautifulSoup
import importlib.resources
import os

SUPPORTED_MATS = [425, 600, 2631, 8237]

def load_cross_section(mat_code, mt_code):
    filename = f"{mat_code}-{mt_code}.html"

    try:
        data_dir = importlib.resources.files("mcmr").parent / "data"
        file_path = os.path.join(data_dir, filename)
    except Exception:
        file_path = ""

    if not os.path.exists(file_path):
        base_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        file_path = os.path.join(base_dir, "data", filename)

    if not os.path.exists(file_path):
        file_path = os.path.join(os.getcwd(), "data", filename)

    if not os.path.exists(file_path):
        raise FileNotFoundError(
            f"File data '{filename}' tidak ditemukan! "
            f"Pastikan folder 'data/' terinstall bersama paket mcmr."
        )

    with open(file_path, "r") as f:
        soup = BeautifulSoup(f, 'html.parser')
    
    lines = soup.find('pre').text.splitlines()
    line = [lines[i].split() for i in range(8, len(lines)-1)]
    df = pd.DataFrame(line[1:], columns=line[0])
    
    E_data = df.iloc[:, 0].astype(float).tolist()
    Sig_data = df.iloc[:, 1].astype(float).tolist()
    
    return E_data, Sig_data


def load_all_materials():
    """Load otomatis seluruh data cross section 4 material."""
    E_total, Sig_total = {}, {}
    E_scatter, Sig_scatter = {}, {}

    for mat in SUPPORTED_MATS:
        e_t, s_t = load_cross_section(mat, 1)
        e_s, s_s = load_cross_section(mat, 2)
        E_total[mat] = e_t
        Sig_total[mat] = s_t
        E_scatter[mat] = e_s
        Sig_scatter[mat] = s_s

    return E_total, Sig_total, E_scatter, Sig_scatter