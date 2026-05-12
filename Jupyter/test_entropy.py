import numpy as np
import matplotlib.pyplot as plt

n_val = 10000
w_val = 0.5
save_filename = f"plot_data_n{n_val}_w{str(w_val).replace('.', '_')}"

arr_T = np.load(f"{save_filename}_T_all_n.npy")
taus = np.load(f"{save_filename}_taus.npy")
beads = np.load(f"{save_filename}_beads.npy")

target_n = 10000
n_idx = target_n - 1
target_bead = 4
bead_idx = np.where(beads == target_bead)[0][0]

E_vals = arr_T[n_idx, bead_idx, :] / target_n

dE = np.diff(E_vals)
tau_mid = (taus[1:] + taus[:-1]) / 2.0

# Since tau is increasing, dE is negative (energy decreases as T decreases / tau increases)
# dS = tau * dE. But we are going from large tau (T=0) down to tau.
# S(tau) = S(tau_max) - \int_{tau}^{tau_max} tau' dE/dtau' dtau'
# So S_i = sum_{j=i}^{N-1} (- tau_mid_j * dE_j)
dS = -tau_mid * dE
S_vals = np.cumsum(dS[::-1])[::-1]

# Let's print some values to check
print(S_vals[:10])
print(S_vals[-10:])
