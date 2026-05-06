import json

with open('BaxterMeth.ipynb', 'r') as f:
    nb = json.load(f)

with open('baxter.jl', 'w') as f:
    for cell in nb['cells']:
        if cell['cell_type'] == 'code':
            f.write("".join(cell['source']) + "\n")
