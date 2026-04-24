import json

with open('JonasPolyConv-SphericalPotential.ipynb', 'r') as f:
    nb = json.load(f)

for cell in nb['cells']:
    if cell['cell_type'] == 'code':
        source = "".join(cell['source'])
        if "plot_runtime" in source:
            new_source = source.replace("d=d", "d=None").replace("d=2", "d=None")
            cell['source'] = new_source.splitlines(keepends=True)

with open('JonasPolyConv-SphericalPotential.ipynb', 'w') as f:
    json.dump(nb, f, indent=1)

