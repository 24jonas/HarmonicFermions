import json

def patch_notebook():
    with open('/home/jonas/Documents/GitHub/HarmonicFermions/Jupyter/BaxterMeth.ipynb', 'r') as f:
        nb = json.load(f)

    # 2. Modify cell 6 (the worker and plot logic)
    source_6 = "".join(nb['cells'][6]['source'])
    
    new_source_6 = """import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from joblib import Parallel, delayed


def energy(n, w, verbose=False):
    \"\"\"
    Computes the exact energy for n fermions in a 2D harmonic trap, 
    scaled by the effective relative frequency w.
    \"\"\"
    # 1. Compute the standard base energy (w=1) using your exact state-filling logic
    cur_lvl = 0
    total_eng = 0
    m = n
    while n > 0:
        for i in range(0, cur_lvl):
            if n > 0:
                total_eng += cur_lvl - 1
                n -= 1
            else:
                break
        cur_lvl += 1
        
    base_energy = total_eng + m   
    
    if verbose:
        print("Highest shell: ", cur_lvl)
        
    # 2. Scale the relative energy component by w
    # For a 2D harmonic oscillator, the 1-particle CM energy is exactly 1.0
    e_1 = 1.0 
    
    true_energy_w = e_1 + w * (base_energy - e_1)
    
    return true_energy_w

# ==========================================
# 1. PA Propagator Functions
# ==========================================
def p_funcs_zeta_1(e):
    return 1.0 + (e**2) / 2.0

def p_funcs_lambda(e):
    return 1.0

def p_funcs_gamma(e):
    return math.sqrt(1.0 + (e**2) / 4.0)

def p_funcs_k1(e):
    return e

# ==========================================
# 2. Factor Calculations
# ==========================================
def factor_calc_T(lambda_val, gamma_val, w, lambda_val_s, gamma_val_s):
    return (lambda_val / gamma_val), (w * lambda_val_s / gamma_val_s)

def factor_calc_H(lambda_val, gamma_val, w, lambda_val_s, gamma_val_s):
    return 0.5 * (gamma_val + 1.0 / gamma_val), (w / 2.0) * (gamma_val_s + 1.0 / gamma_val_s)

# ==========================================
# 3. Top-Level Worker Function (Finite Diff)
# ==========================================
def _compute_estimators_worker_fd(args):
    tau_mid = args[0]
    tau_step = args[1]
    n = args[2]
    N = args[3]
    w = args[4]
    store_all_n = args[5] if len(args) > 5 else False
    
    # Left and right tau for the discrete derivative
    tau_1 = tau_mid - 0.5 * tau_step
    tau_2 = tau_mid + 0.5 * tau_step
    
    # --- A. Evaluate Factors at the midpoint ---
    epsilon = tau_mid / N
    eps_s = w * epsilon
    
    lambda_val = p_funcs_lambda(epsilon)
    gamma_val = p_funcs_gamma(epsilon)
    
    lambda_val_s = p_funcs_lambda(eps_s)
    zeta_1_s = p_funcs_zeta_1(eps_s)
    k1_s = p_funcs_k1(eps_s)
    
    # Guard against domain error in sqrt
    gamma_val_s = math.sqrt(max(0, zeta_1_s**2 - 1.0)) / k1_s
    
    fT_reg, fT_star = factor_calc_T(lambda_val, gamma_val, w, lambda_val_s, gamma_val_s)
    fH_reg, fH_star = factor_calc_H(lambda_val, gamma_val, w, lambda_val_s, gamma_val_s)
    
    # --- B. Compute Energies using YOUR Discrete Derivative ---
    def fd_log_Q(num_particles, b1, b2):
        \"\"\"Finite difference derivative: [log_Q(tau1) - log_Q(tau2)] / step\"\"\"
        lq1 = log_Q(num_particles, b1)
        lq2 = log_Q(num_particles, b2)
        return (lq1 - lq2) / tau_step
        
    def fd_log_Q_all(num_particles, b1, b2):
        lq1 = log_Q_all(num_particles, b1)
        lq2 = log_Q_all(num_particles, b2)
        return (lq1 - lq2) / tau_step

    def get_b_val(tau_val, N_val, w_val):
        eps = w_val * (tau_val / N_val)
        z1 = p_funcs_zeta_1(eps)
        u = math.acosh(z1) if z1 >= 1.0 else 0.0
        return math.exp(-N_val * u)
        
    b_tau1 = get_b_val(tau_1, N, 1.0)
    b_tau2 = get_b_val(tau_2, N, 1.0)
    
    b_s_tau1 = get_b_val(tau_1, N, w)
    b_s_tau2 = get_b_val(tau_2, N, w)

    if store_all_n:
        energy1_T = fd_log_Q(1, b_tau1, b_tau2)
        energystar_T_all = fd_log_Q_all(n, b_s_tau1, b_s_tau2) # shape (n+1,)
        energy1star_T = energystar_T_all[1]
        
        energy_T_all = energy1_T + (energystar_T_all - energy1star_T)
        
        energy1_H = energy1_T * (fH_reg / fT_reg)
        energystar_diff_H_all = (energystar_T_all - energy1star_T) * (fH_star / fT_star)
        energy_H_all = energy1_H + energystar_diff_H_all
        
        # We slice [1:] so index 0 is n=1, shape becomes (n,)
        return energy_T_all[1:], energy_H_all[1:]
    else:
        # Because these use finite difference w.r.t tau, these ARE the Thermodynamic energies!
        energy1_T     = fd_log_Q(1, b_tau1, b_tau2)
        energy1star_T = fd_log_Q(1, b_s_tau1, b_s_tau2)
        energystar_T  = fd_log_Q(n, b_s_tau1, b_s_tau2)
        
        # --- C. Assemble Final Estimators ---
        
        # 1. Total Thermodynamic Estimator
        # No factors needed here, the discrete derivative already contains them!
        energy_T = energy1_T + (energystar_T - energy1star_T)
        
        # 2. Total Hamiltonian Estimator
        # We apply your logic: multiply by inverse thermo factor (1/fT), then by Hamiltonian factor (fH)
        energy1_H = energy1_T * (fH_reg / fT_reg)
        energystar_diff_H = (energystar_T - energy1star_T) * (fH_star / fT_star)
        
        energy_H = energy1_H + energystar_diff_H
        
        return energy_T, energy_H

# ==========================================
# 4. Main Plotting & Saving Function
# ==========================================
def plot_fd_vs_tau_dual(n, N_list, w, tau_start, tau_end, tau_step, save_filename="plot_data", store_all_n=False):
    if tau_step <= 0:
        raise ValueError("tau_step must be positive")
    if tau_end <= tau_start:
        raise ValueError("tau_end must be greater than tau_start")

    taus_left = np.arange(tau_start, tau_end, tau_step)
    taus_left = taus_left[taus_left + tau_step <= tau_end + 1e-15]
    taus_mid = taus_left + 0.5 * tau_step
    
    plt.figure(figsize=(10, 6))
    
    # Initialize dictionaries to hold our results for CSV
    results_T = {"tau_mid": taus_mid}
    results_H = {"tau_mid": taus_mid}
    
    results_T_all_n = []
    results_H_all_n = []
    
    for N in N_list:
        tasks = [(tau, tau_step, n, N, w, store_all_n) for tau in taus_mid]

        # Use joblib to compute estimators in parallel
        results = Parallel(n_jobs=-1)(delayed(_compute_estimators_worker_fd)(task) for task in tasks)
        
        energies_T = np.array([res[0] for res in results])
        energies_H = np.array([res[1] for res in results])
        
        if store_all_n:
            results_T_all_n.append(energies_T)
            results_H_all_n.append(energies_H)
            
            # Plot the highest n (index n-1)
            line, = plt.plot(taus_mid, energies_T[:, -1], label=f"N={N} (Thermo)")
            plt.plot(taus_mid, energies_H[:, -1], linestyle="--", color=line.get_color(), label=f"N={N} (Ham)")
        else:
            results_T[f"N_{N}"] = energies_T
            results_H[f"N_{N}"] = energies_H
            
            # Plot Thermodynamic as solid lines, Hamiltonian as dashed
            line, = plt.plot(taus_mid, energies_T, label=f"N={N} (Thermo)")
            plt.plot(taus_mid, energies_H, linestyle="--", color=line.get_color(), label=f"N={N} (Ham)")

    # Try to plot baseline energy if the `energy(n)` function is defined in your notebook
    try:
        e = energy(n,w)
        print(f"Exact energy for n={n}, w={w}: {e}")
        plt.axhline(e, linestyle=':', color='black', label=f"True Energy (n={n})")
    except NameError:
        pass 

    plt.xlabel(r"$\\tau + \\frac{1}{2}\\,\\Delta\\tau$")
    plt.ylabel(r"Energy Estimators")
    plt.title(f"Thermodynamic (Solid) & Hamiltonian (Dashed) via Finite Difference\\n(n={n}, w={w}, PA Propagator)")
    
    # Put legend outside the plot so it doesn't overlap lines
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    
    # Adjust axes limits as needed based on your n
    try:
        e = energy(n,w)
        plt.ylim(e*0.95, e*1.05) 
    except:
        pass
        
    # plt.xlim(0, 15)
    
    plt.grid(True)
    plt.tight_layout()
    #plt.ylim(4650,5150)
    plt.show()
    
    if store_all_n:
        # Save as npy array with shape (n, num_beads, num_taus)
        arr_T = np.array(results_T_all_n).transpose(2, 0, 1)
        arr_H = np.array(results_H_all_n).transpose(2, 0, 1)
        
        np.save(f"{save_filename}_T_all_n.npy", arr_T)
        np.save(f"{save_filename}_H_all_n.npy", arr_H)
        
        # Save axes to separate files
        np.save(f"{save_filename}_taus.npy", taus_mid)
        np.save(f"{save_filename}_beads.npy", np.array(N_list))
        
        print(f"Data successfully saved to '{save_filename}_T_all_n.npy' (shape {arr_T.shape}) and '{save_filename}_H_all_n.npy' (shape {arr_H.shape})")
        print(f"Axes mapping: n (axis 0), beads (axis 1) -> '{save_filename}_beads.npy', taus (axis 2) -> '{save_filename}_taus.npy'")
        return arr_T, arr_H
    else:
        # Convert to pandas DataFrames and save to CSV
        df_T = pd.DataFrame(results_T)
        df_H = pd.DataFrame(results_H)
        
        df_T.to_csv(f"{save_filename}_T.csv", index=False)
        df_H.to_csv(f"{save_filename}_H.csv", index=False)
        
        print(f"Data successfully saved to '{save_filename}_T.csv' and '{save_filename}_H.csv'")
        
        return df_T, df_H
"""

    # We need to split new_source_6 by lines and write as array to match notebook format
    nb['cells'][6]['source'] = [line + '\n' for line in new_source_6.split('\n')]
    # fix the last newline
    nb['cells'][6]['source'][-1] = nb['cells'][6]['source'][-1].rstrip('\n')

    with open('/home/jonas/Documents/GitHub/HarmonicFermions/Jupyter/BaxterMeth.ipynb', 'w') as f:
        json.dump(nb, f, indent=1)

if __name__ == "__main__":
    patch_notebook()
