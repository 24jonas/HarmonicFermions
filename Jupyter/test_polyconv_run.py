import json

notebook_path = "/home/jonas/Documents/GitHub/HarmonicFermions/Jupyter/JonasPolyConv-opt.ipynb"

with open(notebook_path, "r", encoding="utf-8") as f:
    nb = json.load(f)

# Extract cell 1 and cell 3
code = "".join(nb['cells'][2]['source']) + "\n" + "".join(nb['cells'][4]['source'])
with open("temp_polyconv.py", "w") as f:
    f.write(code)
    f.write("\nprint('Testing fermion_logZ_numeric...')\n")
    f.write("res = fermion_logZ_numeric(tau=10.0, N=200, n=7, d=2)\n")
    f.write("print('Result logZ_n:', res)\n")

