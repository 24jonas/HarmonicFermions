import json

def append_extraction_cell():
    with open('/home/jonas/Documents/GitHub/HarmonicFermions/Jupyter/BaxterMeth.ipynb', 'r') as f:
        nb = json.load(f)

    # Define the new cell source
    new_cell_source = """import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

# Load the saved 3D arrays and axes
# Change save_filename if you used a different name in plot_fd_vs_tau_dual
save_filename = "plot_data"  

try:
    arr_T = np.load(f"{save_filename}_T_all_n.npy")
    arr_H = np.load(f"{save_filename}_H_all_n.npy")
    taus = np.load(f"{save_filename}_taus.npy")
    beads = np.load(f"{save_filename}_beads.npy")
    
    print(f"Loaded arrays with shapes:")
    print(f"T_all_n: {arr_T.shape}")
    print(f"H_all_n: {arr_H.shape}")
    print(f"taus: {taus.shape}")
    print(f"beads: {beads.shape}")
    
    # ==========================================
    # Define the target n and bead you want to extract
    # ==========================================
    target_n = 2   # Example: n=2
    target_bead = 4  # Example: N=4
    
    # 1. Find the correct indices
    # Remember that array index 0 corresponds to n=1, so index = target_n - 1
    n_idx = target_n - 1
    
    # Find the index in the beads array that matches target_bead
    bead_idx = np.where(beads == target_bead)[0][0]
    
    print(f"\\nExtracting data for n={target_n} (index {n_idx}) and bead={target_bead} (index {bead_idx})")
    
    # 2. Extract the 1D arrays for tau vs energy
    energy_T = arr_T[n_idx, bead_idx, :]
    energy_H = arr_H[n_idx, bead_idx, :]
    
    # 3. Save this specific 1D slice to a CSV file
    df_extracted = pd.DataFrame({
        "tau_mid": taus,
        "Energy_Thermodynamic": energy_T,
        "Energy_Hamiltonian": energy_H
    })
    export_filename = f"extracted_n{target_n}_bead{target_bead}.csv"
    df_extracted.to_csv(export_filename, index=False)
    print(f"Successfully exported to {export_filename}\\n")
    
    # 4. Plot to verify
    plt.figure(figsize=(8, 5))
    plt.plot(taus, energy_T, label=f"N={target_bead} (Thermo)", linestyle="-")
    plt.plot(taus, energy_H, label=f"N={target_bead} (Ham)", linestyle="--")
    
    # Try to plot baseline energy if the `energy(n)` function is defined
    try:
        e = energy(target_n, w=1.0)
        plt.axhline(e, linestyle=':', color='black', label=f"True Energy (n={target_n})")
    except NameError:
        pass 
    
    plt.xlabel(r"$\\tau + \\frac{1}{2}\\,\\Delta\\tau$")
    plt.ylabel(r"Energy")
    plt.title(f"Extracted Data: n={target_n}, Bead N={target_bead}")
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.show()

except FileNotFoundError:
    print(f"Error: Could not find the files starting with '{save_filename}'.")
    print("Please make sure you have run plot_fd_vs_tau_dual with store_all_n=True first!")
"""
    
    # Format the source properly for Jupyter (list of strings with newlines)
    source_lines = [line + '\\n' for line in new_cell_source.split('\\n')]
    if source_lines:
        source_lines[-1] = source_lines[-1].rstrip('\\n')
        
    new_cell = {
        "cell_type": "code",
        "execution_count": None,
        "metadata": {},
        "outputs": [],
        "source": source_lines
    }
    
    # Append the new cell to the notebook
    nb['cells'].append(new_cell)

    with open('/home/jonas/Documents/GitHub/HarmonicFermions/Jupyter/BaxterMeth.ipynb', 'w') as f:
        json.dump(nb, f, indent=1)

if __name__ == "__main__":
    append_extraction_cell()
