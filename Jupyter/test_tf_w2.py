import numpy as np
from scipy.integrate import quad
from scipy.optimize import root_scalar

n_target = 10000
w_val = 0.5
save_filename = f"plot_data_n{n_target}_w{str(w_val).replace('.', '_')}"

arr_T = np.load(f"{save_filename}_T_all_n.npy")
taus = np.load(f"{save_filename}_taus.npy")
beads = np.load(f"{save_filename}_beads.npy")

target_bead = 4
target_tau_idx = 50
bead_idx = np.where(beads == target_bead)[0][0]
tau = taus[target_tau_idx]
T_fixed = 1.0 / tau

E_T_n = arr_T[n_target-1, bead_idx, target_tau_idx] / n_target

def N_int(mu, T):
    def integrand(e):
        x = (e - mu)/T
        if x > 100: return 0.0
        return (e / (w_val**2) + 1.0 / w_val) / (np.exp(x) + 1.0)
    return quad(integrand, 0, max(0, mu) + 40*T)[0]

def E_int(mu, T):
    def integrand(e):
        x = (e - mu)/T
        if x > 100: return 0.0
        return (e**2 / (w_val**2) + e / w_val) / (np.exp(x) + 1.0)
    return quad(integrand, 0, max(0, mu) + 40*T)[0]

mu0 = np.sqrt(2 * n_target) * w_val

def obj1(mu): return N_int(mu, T_fixed) - n_target
res1 = root_scalar(obj1, bracket=[-20*T_fixed, mu0 + 20*T_fixed], method='brentq')
mu_TF = res1.root
E_TF_n = E_int(mu_TF, T_fixed) / n_target

print(f"Simulated E/n: {E_T_n}")
print(f"TF E/n (with 1/w term): {E_TF_n}")

