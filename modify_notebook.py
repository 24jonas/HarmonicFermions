import json

with open('/home/jonas/Documents/GitHub/HarmonicFermions/Jupyter/BaxterMeth.ipynb', 'r') as f:
    nb = json.load(f)

for i, cell in enumerate(nb['cells']):
    if cell['cell_type'] == 'code':
        source = "".join(cell['source'])
        print(f"Cell {i} length: {len(source)}")
        if "def log_Q(n, b):" in source:
            print(f"-> Found log_Q in cell {i}")
        if "def _compute_estimators_worker_fd(args):" in source:
            print(f"-> Found _compute_estimators_worker_fd in cell {i}")
