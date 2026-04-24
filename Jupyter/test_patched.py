import json

with open('JonasPolyConv-Copy3.ipynb', 'r') as f:
    nb = json.load(f)

code = ""
for cell in nb['cells']:
    if cell['cell_type'] == 'code':
        source = "".join(cell['source'])
        code += source + "\n"

# write to a python file
with open('extracted3.py', 'w') as f:
    f.write(code)

