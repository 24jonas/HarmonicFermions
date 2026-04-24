import json

with open('JonasPolyConv-SphericalPotential.ipynb', 'r') as f:
    nb = json.load(f)

for cell in nb['cells']:
    if cell['cell_type'] == 'code':
        source = "".join(cell['source'])
        if "results = plot_logZ_and_fd_multiN(" in source:
            new_source = source.replace("N_list=N,", "N_list=[N],")
            cell['source'] = new_source.splitlines(keepends=True)

with open('JonasPolyConv-SphericalPotential.ipynb', 'w') as f:
    json.dump(nb, f, indent=1)

