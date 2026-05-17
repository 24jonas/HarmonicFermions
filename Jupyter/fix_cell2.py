import json

path = "/home/jonas/Documents/GitHub/HarmonicFermions/Jupyter/heat_capacity.ipynb"
with open(path, "r") as f:
    nb = json.load(f)

new_source = """# --- Comparison with Thomas-Fermi (Semiclassical) Limit vs n ---
from scipy.integrate import quad
from scipy.optimize import root_scalar
import time
import os
import numpy as np
import matplotlib.pyplot as plt

print("Computing Thomas-Fermi limit vs n... (this might take a few seconds)")

target_bead = 4
target_tau_idx = 100  # Choose a specific tau index for the plots

bead_idx = np.where(beads == target_bead)[0][0]
tau = taus[target_tau_idx]
T_fixed = 1.0 / tau

# 1. Extract Energy per fermion vs n
n_vals_all = np.arange(1, n_val + 1)
E_T_n = arr_T[:, bead_idx, target_tau_idx] / n_vals_all
E_H_n = arr_H[:, bead_idx, target_tau_idx] / n_vals_all

target_bead_high = 256
has_high_beads = False
if 256 in beads:
    exact_bead_idx = np.where(beads == 256)[0][0]
    n_val_high = arr_H.shape[0]
    E_H_n_high = arr_H[:, exact_bead_idx, target_tau_idx] / np.arange(1, n_val_high + 1)
    has_high_beads = True
else:
    print("Could not find high bead data in loaded arrays.")

# 2. Compute Heat Capacity vs n
# We differentiate along the tau axis (-1) for all n at once
dE_T_tau = np.diff(arr_T[:, bead_idx, :], axis=-1)
dE_H_tau = np.diff(arr_H[:, bead_idx, :], axis=-1)
tau_mid = (taus[1:] + taus[:-1]) / 2.0

Cv_T_all = -(tau_mid**2) * (dE_T_tau / np.diff(taus))
Cv_H_all = -(tau_mid**2) * (dE_H_tau / np.diff(taus))

if has_high_beads:
    dE_H_tau_high = np.diff(arr_H[:, exact_bead_idx, :], axis=-1)
    Cv_H_all_high = -(tau_mid**2) * (dE_H_tau_high / np.diff(taus))

# Note: target_tau_idx might need adjusting because tau_mid is 1 element shorter.
Cv_T_n = Cv_T_all[:, target_tau_idx] / n_vals_all
Cv_H_n = Cv_H_all[:, target_tau_idx] / n_vals_all
if has_high_beads:
    Cv_H_n_high = Cv_H_all_high[:, target_tau_idx] / np.arange(1, n_val_high + 1)

# 3. Compute Chemical Potential vs n
mu_T_all = np.diff(arr_T[:, bead_idx, target_tau_idx])
mu_H_all = np.diff(arr_H[:, bead_idx, target_tau_idx])
n_vals_mu = np.arange(2, n_val + 1)
if has_high_beads:
    mu_H_all_high = np.diff(arr_H[:, exact_bead_idx, target_tau_idx])

# --- Thomas-Fermi Calculations vs n ---
# Ensure we use the correct frequency for the TF integrals
# w_val is already set in the first cell, so we use it directly.
def N_int(mu, T):
    def integrand(e):
        x = (e - mu)/T
        if x > 100: return 0.0
        return (e / w_val**2) / (np.exp(x) + 1.0)
    return quad(integrand, 0, max(0, mu) + 40*T)[0]

def E_int(mu, T):
    def integrand(e):
        x = (e - mu)/T
        if x > 100: return 0.0
        return (e**2 / w_val**2) / (np.exp(x) + 1.0)
    return quad(integrand, 0, max(0, mu) + 40*T)[0]

# Calculate TF limit for a sparse set of n values to make it fast and plot as a smooth line
n_vals_TF = np.linspace(10, n_val, 50)
E_TF_n = np.zeros_like(n_vals_TF)
Cv_TF_n = np.zeros_like(n_vals_TF)
mu_TF = np.zeros_like(n_vals_TF)

# To compute TF Cv, we need E_TF at T and T+dT
dT = 0.01
T2 = T_fixed + dT

for i, n_target in enumerate(n_vals_TF):
    mu0 = np.sqrt(2 * n_target) * w_val
    
    # At T_fixed
    def obj1(mu): return N_int(mu, T_fixed) - n_target
    res1 = root_scalar(obj1, bracket=[-20*T_fixed, mu0 + 20*T_fixed], method='brentq')
    mu_TF[i] = res1.root
    E1 = E_int(mu_TF[i], T_fixed)
    E_TF_n[i] = E1 / n_target
    
    # At T2 for Cv
    def obj2(mu): return N_int(mu, T2) - n_target
    res2 = root_scalar(obj2, bracket=[-20*T2, mu0 + 20*T2], method='brentq')
    E2 = E_int(res2.root, T2)
    
    Cv_TF_n[i] = (E2 - E1) / dT / n_target

# --- Plots ---

ntop = 1500

plt.figure(figsize=(10, 6))
plt.xlim(0,ntop)
plt.ylim(0,20)
plt.plot(n_vals_all, E_T_n, label=f"Thermodynamic (N={target_bead})", linestyle="-")
plt.plot(n_vals_all, E_H_n, label=f"Hamiltonian (N={target_bead})", linestyle="--")
if has_high_beads:
    plt.plot(np.arange(1, n_val_high + 1), E_H_n_high, label=f"Exact N={target_bead_high} (Ham)", linestyle="--", color="black", linewidth=2)
plt.plot(n_vals_TF, E_TF_n, label="Thomas-Fermi", linestyle=":", color="black", linewidth=2)
plt.xlabel("Number of Fermions (n)")
plt.ylabel(r"Energy per Fermion ($E / n$)")
plt.title(f"Energy vs n ($T$={T_fixed:.3f}, w={w_val})")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()

plt.figure(figsize=(10, 6))
plt.xlim(0,ntop)
plt.ylim(-20,20)
plt.plot(n_vals_all, Cv_T_n, label=f"Thermodynamic (N={target_bead})", linestyle="-")
plt.plot(n_vals_all, Cv_H_n, label=f"Hamiltonian (N={target_bead})", linestyle="--")
if has_high_beads:
    plt.plot(np.arange(1, n_val_high + 1), Cv_H_n_high, label=f"Exact N={target_bead_high} (Ham)", linestyle="--", color="black", linewidth=2)
plt.plot(n_vals_TF, Cv_TF_n, label="Thomas-Fermi", linestyle=":", color="black", linewidth=2)
plt.xlabel("Number of Fermions (n)")
plt.ylabel(r"Heat Capacity per Fermion ($C_v / n$)")
plt.title(f"Heat Capacity vs n ($T_{{mid}}$={1.0/tau_mid[target_tau_idx]:.3f}, w={w_val})")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()

plt.figure(figsize=(10, 6))
plt.xlim(0,ntop)
plt.ylim(0,40)
plt.plot(n_vals_mu, mu_T_all, label=f"Thermodynamic (N={target_bead})", linestyle="-")
plt.plot(n_vals_mu, mu_H_all, label=f"Hamiltonian (N={target_bead})", linestyle="--")
if has_high_beads:
    plt.plot(np.arange(2, n_val_high + 1), mu_H_all_high, label=f"Exact N={target_bead_high} (Ham)", linestyle="--", color="black", linewidth=2)
plt.plot(n_vals_TF, mu_TF, label="Thomas-Fermi", linestyle=":", color="black", linewidth=2)
plt.xlabel("Number of Fermions (n)")
plt.ylabel(r"Chemical Potential ($\mu$)")
plt.title(f"Chemical Potential vs n ($T$={T_fixed:.3f}, w={w_val})")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()
"""

new_lines = [line + "\n" for line in new_source.split("\n")]
# Remove the last newline from the very last empty line to avoid adding extra blank line
if new_lines and new_lines[-1] == "\n":
    new_lines.pop()

nb["cells"][2]["source"] = new_lines

with open(path, "w") as f:
    json.dump(nb, f, indent=1)

