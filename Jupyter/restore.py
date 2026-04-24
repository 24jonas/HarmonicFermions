import json

with open('JonasPolyConv.ipynb', 'r') as f:
    nb_orig = json.load(f)

for cell in nb_orig['cells']:
    if cell['cell_type'] == 'code':
        source = "".join(cell['source'])
        if "def plot_logZ_and_fd_multiN" in source:
            def_source = cell['source']
            break

with open('JonasPolyConv-Copy3.ipynb', 'r') as f:
    nb = json.load(f)

# Find the cell that has "results = plot_logZ_and_fd_multiN(" at the top (which used to be the definition)
for i, cell in enumerate(nb['cells']):
    if cell['cell_type'] == 'code':
        source = "".join(cell['source'])
        if source.startswith("n = 200\nresults = plot_logZ_and_fd_multiN"):
            # this is the first overwritten one. We'll restore it to def_source.
            cell['source'] = def_source
            break

with open('JonasPolyConv-Copy3.ipynb', 'w') as f:
    json.dump(nb, f, indent=1)

