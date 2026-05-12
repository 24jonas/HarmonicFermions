import numpy as np

n_target = 10000
w_val = 0.5
save_filename = f"plot_data_n{n_target}_w{str(w_val).replace('.', '_')}"

arr_T = np.load(f"{save_filename}_T_all_n.npy")
beads = np.load(f"{save_filename}_beads.npy")
target_tau_idx = 50

for b_idx, b in enumerate(beads):
    E_T_n = arr_T[n_target-1, b_idx, target_tau_idx] / n_target
    print(f"Bead={b}: Simulated E/n = {E_T_n}")
