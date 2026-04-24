import json

with open('JonasPolyConv-Copy3.ipynb', 'r') as f:
    nb = json.load(f)

for cell in nb['cells']:
    if cell['cell_type'] == 'code':
        source = "".join(cell['source'])
        if "results = plot_logZ_and_fd_multiN(" in source and "tau_end=50," in source:
            cell['source'] = [line.replace("tau_end=50,", "tau_end=5000,") for line in cell['source']]

with open('JonasPolyConv-Copy3.ipynb', 'w') as f:
    json.dump(nb, f, indent=1)

