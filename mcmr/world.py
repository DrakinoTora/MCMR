import xml.etree.ElementTree as ET

__all__ = ["export_world", "load_world"]


def export_world(filename, x_world, y_world, x_grid, y_grid, material_matrix, sources,
                  bc_top="vacuum", bc_bot="vacuum", bc_left="vacuum", bc_right="vacuum"):
    """Simpan definisi dunia simulasi (geometri, grid, material, sumber, boundary) ke XML.

    File ini TERPISAH dari XML hasil simulasi (mcmr_results.xml) -- 1 file = 1 informasi.
    Cukup di-export sekali, lalu dipakai berkali-kali lewat mcmr.run_from_world() dan
    ResultsPlotter tanpa perlu mengetik ulang parameter grid.

    material_matrix, sources : list 2D dengan konvensi [row][col] sama seperti load_and_run:
                         row=0 paling ATAS (y tertinggi), col=0 paling KIRI (x=0).

    Return: filename yang ditulis.
    """
    ny = len(material_matrix)
    nx = len(material_matrix[0]) if ny else 0
    if any(len(row) != nx for row in material_matrix):
        raise ValueError("setiap baris material_matrix harus punya jumlah kolom yang sama")
    if len(sources) != ny or any(len(row) != nx for row in sources):
        raise ValueError("sources harus punya shape yang sama persis dengan material_matrix")

    root = ET.Element("mcmr_world")

    dims = ET.SubElement(root, "dimensions")
    dims.set("x_world", str(x_world))
    dims.set("y_world", str(y_world))

    grid_el = ET.SubElement(root, "grid")
    ET.SubElement(grid_el, "x_grid").text = ",".join(map(str, x_grid))
    ET.SubElement(grid_el, "y_grid").text = ",".join(map(str, y_grid))

    bc_el = ET.SubElement(root, "boundary")
    bc_el.set("top", bc_top)
    bc_el.set("bottom", bc_bot)
    bc_el.set("left", bc_left)
    bc_el.set("right", bc_right)

    mats_el = ET.SubElement(root, "materials")
    for i, row in enumerate(material_matrix):
        row_el = ET.SubElement(mats_el, "row")
        row_el.set("index", str(i))
        row_el.text = ",".join(row)

    src_el = ET.SubElement(root, "sources")
    for i, row in enumerate(sources):
        row_el = ET.SubElement(src_el, "row")
        row_el.set("index", str(i))
        row_el.text = ",".join(map(str, row))

    tree = ET.ElementTree(root)
    ET.indent(tree, space="  ")
    tree.write(filename, xml_declaration=True, encoding="UTF-8")
    return filename


def load_world(filename):
    """Baca file world XML hasil export_world(), kembalikan dict siap-pakai.

    Keys: x_world, y_world, x_grid, y_grid, material_matrix, sources,
          bc_top, bc_bot, bc_left, bc_right
    """
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

    return {
        "x_world": x_world, "y_world": y_world,
        "x_grid": x_grid, "y_grid": y_grid,
        "material_matrix": material_matrix, "sources": sources,
        "bc_top": bc_top, "bc_bot": bc_bot, "bc_left": bc_left, "bc_right": bc_right,
    }