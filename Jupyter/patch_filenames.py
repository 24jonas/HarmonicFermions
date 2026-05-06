import json

def patch_baxter():
    with open("BaxterMeth-opt.ipynb", "r") as f:
        nb = json.load(f)
    
    cell2 = "".join(nb["cells"][2]["source"])
    cell2 = cell2.replace('save_filename="plot_data"', 'save_filename=f"plot_data_n{n}_w{str(w).replace(\'.\', \'_\')}"')
    nb["cells"][2]["source"] = [line + ("\n" if not line.endswith("\n") else "") for line in cell2.splitlines()]
    
    cell3 = "".join(nb["cells"][3]["source"])
    cell3 = cell3.replace('save_filename = "plot_data"', 'n_val = 50050\nw_val = 1.0\nsave_filename = f"plot_data_n{n_val}_w{str(w_val).replace(\'.\', \'_\')}"')
    nb["cells"][3]["source"] = [line + ("\n" if not line.endswith("\n") else "") for line in cell3.splitlines()]
    
    with open("BaxterMeth-opt.ipynb", "w") as f:
        json.dump(nb, f, indent=1)
        
def patch_jonas():
    with open("JonasPolyConv-opt.ipynb", "r") as f:
        nb = json.load(f)
        
    cell2 = "".join(nb["cells"][2]["source"])
    if "n_val" not in cell2:
        cell2 = cell2.replace('n=300,', 'n=n_val,')
        cell2 = cell2.replace('w=1.0,', 'w=w_val,')
        cell2 = "n_val = 300\nw_val = 1.0\n\n" + cell2
    cell2 = cell2.replace('save_filename="plot_data"', 'save_filename=f"plot_data_n{n_val}_w{str(w_val).replace(\'.\', \'_\')}"')
    nb["cells"][2]["source"] = [line + ("\n" if not line.endswith("\n") else "") for line in cell2.splitlines()]
    
    cell3 = "".join(nb["cells"][3]["source"])
    cell3 = cell3.replace('save_filename = "plot_data"', 'n_val = 300\nw_val = 1.0\nsave_filename = f"plot_data_n{n_val}_w{str(w_val).replace(\'.\', \'_\')}"')
    nb["cells"][3]["source"] = [line + ("\n" if not line.endswith("\n") else "") for line in cell3.splitlines()]
    
    with open("JonasPolyConv-opt.ipynb", "w") as f:
        json.dump(nb, f, indent=1)

patch_baxter()
patch_jonas()
print("Patched filenames successfully.")
