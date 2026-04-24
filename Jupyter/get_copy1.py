import json

with open('JonasPolyConv-Copy1.ipynb', 'r') as f:
    nb = json.load(f)

for cell in nb['cells']:
    if cell['cell_type'] == 'code':
        source = "".join(cell['source'])
        if "def " in source and "val" in source:
            print(source)

