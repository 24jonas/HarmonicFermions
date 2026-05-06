import json

with open("BaxterMeth-opt.ipynb", "r") as f:
    nb = json.load(f)

for i, cell in enumerate(nb["cells"]):
    if cell["cell_type"] == "code":
        src = "".join(cell["source"])
        if "results =" in src or "save_filename" in src:
            print(f"--- Cell {i} ---")
            print(src[:250])
            print("...")
