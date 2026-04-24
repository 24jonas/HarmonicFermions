import json

with open('JonasPolyConv-Copy3.ipynb', 'r') as f:
    nb = json.load(f)

for cell in nb['cells']:
    if cell['cell_type'] == 'code':
        source = "".join(cell['source'])
        if "val =" in source:
            print("FOUND VAL CELL:")
            print(source)
