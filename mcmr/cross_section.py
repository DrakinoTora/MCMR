import pandas as pd
from bs4 import BeautifulSoup
import importlib.resources
import os

def load_cross_section(material_code, mt_code):
    filename = f"{material_code}-{mt_code}.html"
    
    # Cari path file HTML dari folder data
    try:
        data_dir = importlib.resources.files("mcmr").parent / "data"
        file_path = os.path.join(data_dir, filename)
    except Exception:
        base_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        file_path = os.path.join(base_path, "data", filename)

    with open(file_path, "r") as f:
        soup = BeautifulSoup(f, 'html.parser')
    
    lines = soup.find('pre').text.splitlines()
    line = [lines[i].split() for i in range(8, len(lines)-1)]
    df = pd.DataFrame(line[1:], columns=line[0])
    
    E_data = df.iloc[:, 0].astype(float).tolist()
    Sig_data = df.iloc[:, 1].astype(float).tolist()
    
    return E_data, Sig_data