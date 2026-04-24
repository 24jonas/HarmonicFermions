import json

with open('JonasPolyConv-Copy3.ipynb', 'r') as f:
    nb = json.load(f)

for cell in nb['cells']:
    if cell['cell_type'] == 'code':
        source = "".join(cell['source'])
        
        if "val =" in source and "fermion_logZ_numeric" in source:
            new_source = """n = 100

L_box = 2.0
def Ek_box(k):
    return (math.pi**2 / (2 * L_box**2)) * k

val = (
    fermion_logZ_numeric(tau=2.0,  N=200000, n=n, d=2, k_start=1, Ek_fn=Ek_box)
    - fermion_logZ_numeric(tau=2.5,  N=200000, n=n, d=2, k_start=1, Ek_fn=Ek_box)
) / 0.5

print(val)
"""
            # split lines and keep formatting
            lines = [line + "\n" for line in new_source.split('\n')]
            # remove last newline from last line to match standard notebook format
            if lines: lines[-1] = lines[-1][:-1]
            cell['source'] = lines
            cell['outputs'] = [] # clear outputs
            
        elif "plot_logZ_and_fd_multiN" in source and "show_fd=True" in source:
            new_source = """n = 200
results = plot_logZ_and_fd_multiN(
    tau_start=0.1,
    tau_end=50,
    tau_step=2.5,
    
    n=n,
    d=2,
    N_list=[1000,10000,100000],
    show_logZ=False,
    show_fd=True,
    Ek_fn=Ek_box
)
"""
            # NOTE: changed tau_end to 50 instead of 5000 so it can be evaluated in reasonable time later
            lines = [line + "\n" for line in new_source.split('\n')]
            if lines: lines[-1] = lines[-1][:-1]
            cell['source'] = lines
            cell['outputs'] = [] # clear outputs

with open('JonasPolyConv-Copy3.ipynb', 'w') as f:
    json.dump(nb, f, indent=1)

