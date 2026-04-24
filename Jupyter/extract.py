import json

with open('JonasPolyConv-Copy3.ipynb', 'r') as f:
    nb = json.load(f)

with open('extracted2.py', 'w') as f:
    for cell in nb['cells']:
        if cell['cell_type'] == 'code':
            f.write("".join(cell['source']) + "\n")
