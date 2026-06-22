import math
import numpy as np
import pandas as pd
import numba
from baxter_method import log_Q_all
from propagator import p_funcs_zeta_1, p_funcs_lambda, p_funcs_gamma, p_funcs_k1, factor_calc_T, factor_calc_H

@numba.njit
def logsumexp2(x, y):
    if x == -np.inf:
        return y
    if y == -np.inf:
        return x
    m = x if x > y else y
    return m + math.log(math.exp(x - m) + math.exp(y - m))

@numba.njit
def log_Z_B_all(n, b):
    log_Z = np.full(n + 1, -np.inf, dtype=np.float64)
    log_Z[0] = 0.0
    
    logb = math.log(b)
    log_Z1_arr = np.full(n + 1, -np.inf, dtype=np.float64)
    for k in range(1, n + 1):
        bk = math.exp(k * logb)
        log_Z1_arr[k] = k * logb - 2.0 * math.log1p(-bk)
        
    for i in range(1, n + 1):
        s = -np.inf
        for k in range(1, i + 1):
            term = log_Z1_arr[k] + log_Z[i - k]
            s = logsumexp2(s, term)
        log_Z[i] = s - math.log(i)
        
    return log_Z

def get_b_val(tau_val, N_val, w_val):
    if N_val == 0:
        return math.exp(-w_val * abs(tau_val))
    eps = w_val * (tau_val / N_val)
    z1 = p_funcs_zeta_1(eps)
    u = math.acosh(z1) if z1 >= 1.0 else 0.0
    return math.exp(-N_val * u)

def get_log_Z_ref(n_max, b_tau, predict_with):
    if predict_with == 'Z_B':
        return log_Z_B_all(n_max, b_tau)
    elif predict_with == 'Z_D':
        Z1 = b_tau / (1.0 - b_tau)**2
        log_Z1 = math.log(Z1)
        lz_ref = np.zeros(n_max + 1)
        for n_val in range(n_max + 1):
            lz_ref[n_val] = n_val * log_Z1 - math.lgamma(n_val + 1)
        return lz_ref
    else:
        raise ValueError("predict_with must be 'Z_B' or 'Z_D'")

def compute_energies(n_max, tau, N, w, predict_with):
    step = 1e-5
    tau1 = tau - step / 2.0
    tau2 = tau + step / 2.0
    
    b_tau1 = get_b_val(tau1, N, w)
    b_tau2 = get_b_val(tau2, N, w)
    
    def fd_diff(arr1, arr2):
        return (arr1 - arr2) / step
        
    lq1_F = log_Q_all(n_max, b_tau1)
    lq2_F = log_Q_all(n_max, b_tau2)
    energystar_T_all_F = fd_diff(lq1_F, lq2_F)
    
    energy1star_T = energystar_T_all_F[1]
    
    lz1_ref = get_log_Z_ref(n_max, b_tau1, predict_with)
    lz2_ref = get_log_Z_ref(n_max, b_tau2, predict_with)
    energystar_T_all_ref = fd_diff(lz1_ref, lz2_ref)
    
    b_true_tau1 = math.exp(-w * tau1)
    b_true_tau2 = math.exp(-w * tau2)
    
    lq1_1 = log_Q_all(1, b_true_tau1)
    lq2_1 = log_Q_all(1, b_true_tau2)
    energy1_T = fd_diff(lq1_1, lq2_1)[1]
    
    eps = w * tau / N if N > 0 else 0.0
    z1 = p_funcs_zeta_1(eps)
    lam = p_funcs_lambda(eps)
    gam = p_funcs_gamma(eps)
    
    eps_s = w * tau if N > 0 else 0.0
    z1_s = p_funcs_zeta_1(eps_s)
    lam_s = p_funcs_lambda(eps_s)
    gam_s = p_funcs_gamma(eps_s)
    
    fT_reg, fT_star = factor_calc_T(lam, gam, w, lam_s, gam_s)
    fH_reg, fH_star = factor_calc_H(lam, gam, w, lam_s, gam_s)
    
    energy_T_F = energy1_T + (energystar_T_all_F - energy1star_T)
    energy_T_ref = energy1_T + (energystar_T_all_ref - energy1star_T)
    
    energy1_H = energy1_T * (fH_reg / fT_reg)
    energy_H_F = energy1_H + (energystar_T_all_F - energy1star_T) * (fH_star / fT_star)
    energy_H_ref = energy1_H + (energystar_T_all_ref - energy1star_T) * (fH_star / fT_star)
    
    energy_T_F[0] = 0.0
    energy_T_ref[0] = 0.0
    energy_H_F[0] = 0.0
    energy_H_ref[0] = 0.0
    
    return energy_T_F, energy_T_ref, energy_H_F, energy_H_ref

d = 2
tau = 5.0
N = 0
w = 1.0
n_max = 14

predict_with = 'Z_B'  # Options: 'Z_B', 'Z_D'

energy_T_F, energy_T_ref, energy_H_F, energy_H_ref = compute_energies(n_max, tau, N, w, predict_with)

b_tau = get_b_val(tau, N, w)
lq_F = log_Q_all(n_max, b_tau)
lz_ref = get_log_Z_ref(n_max, b_tau, predict_with)

# Correct sign calculation using the free energy / partition functions
# In statistical mechanics, <sgn> = Z_F / Z_ref = exp(ln Z_F - ln Z_ref)
# using E_F - E_ref omits the entropy difference (since F = E - TS)
sgn_true = np.exp(lq_F - lz_ref)

# The old calculation based purely on energy
sgn_E_T = np.exp(-tau * (energy_T_F - energy_T_ref))
sgn_E_H = np.exp(-tau * (energy_H_F - energy_H_ref))

n_arr = np.arange(n_max + 1)
df = pd.DataFrame({
    'n': n_arr,
    'E_f_T': energy_T_F,
    'E_ref_T': energy_T_ref,
    'sgn_true': sgn_true,
    'sgn_E_T': sgn_E_T,
    'E_f_H': energy_H_F,
    'E_ref_H': energy_H_ref,
    'sgn_E_H': sgn_E_H
})

print(df.head(10))
csv_filename = f'sgn_vs_n_data_tau{tau}_N{N}_w{w}.csv'
df.to_csv(csv_filename, index=False)
