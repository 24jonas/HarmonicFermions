import numpy as np

n_target = 10000
w_val = 0.5
save_filename = f"plot_data_n{n_target}_w{str(w_val).replace('.', '_')}"

arr_T = np.load(f"{save_filename}_T_all_n.npy")
arr_H = np.load(f"{save_filename}_H_all_n.npy")
taus = np.load(f"{save_filename}_taus.npy")
beads = np.load(f"{save_filename}_beads.npy")

target_bead = 4
bead_idx = np.where(beads == target_bead)[0][0]

E_T = arr_T[n_target-1, bead_idx, :] / n_target
E_H = arr_H[n_target-1, bead_idx, :] / n_target

print("tau (last 5):", taus[-5:])
print("E_T (last 5):", E_T[-5:])
print("E_H (last 5):", E_H[-5:])
print("dE_T (last 5):", np.diff(E_T)[-5:])
print("dE_H (last 5):", np.diff(E_H)[-5:])
