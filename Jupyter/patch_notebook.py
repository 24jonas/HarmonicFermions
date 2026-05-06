import json

with open("JonasPolyConv-opt.ipynb", "r") as f:
    nb = json.load(f)

def read_cell(filename):
    with open(filename, "r") as f:
        lines = f.readlines()
    return lines

nb["cells"][0]["source"] = read_cell("/home/jonas/.gemini/antigravity/brain/46091432-a5cd-4204-95ec-50eb50161a28/scratch/new_cell0.py")
nb["cells"][1]["source"] = read_cell("/home/jonas/.gemini/antigravity/brain/46091432-a5cd-4204-95ec-50eb50161a28/scratch/new_cell1.py")
nb["cells"][2]["source"] = read_cell("/home/jonas/.gemini/antigravity/brain/46091432-a5cd-4204-95ec-50eb50161a28/scratch/new_cell2.py")
nb["cells"][3]["source"] = read_cell("/home/jonas/.gemini/antigravity/brain/46091432-a5cd-4204-95ec-50eb50161a28/scratch/new_cell3.py")

with open("JonasPolyConv-opt.ipynb", "w") as f:
    json.dump(nb, f, indent=1)

print("Notebook patched successfully.")
