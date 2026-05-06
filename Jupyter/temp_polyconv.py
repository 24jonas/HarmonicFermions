import math
import numpy as np
import numba
import sqlite3
import os
import time

def init_polyconv_db(db_name='recursions_polyconv.db'):
    try:
        conn = sqlite3.connect(db_name, timeout=60.0)
        c = conn.cursor()
        c.execute('''CREATE TABLE IF NOT EXISTS polyconv_states
                     (n INTEGER, d INTEGER, b REAL, k_start INTEGER, logZ BLOB, prev_logZn REAL, small_count INTEGER, cumulative_capacity REAL,
                     PRIMARY KEY (n, d, b))''')
        conn.commit()
        conn.close()
    except Exception:
        pass

def save_polyconv_state(n, d, b, k_start, logZ, prev_logZn, small_count, cumulative_capacity, db_name='recursions_polyconv.db'):
    for _ in range(5):
        try:
            conn = sqlite3.connect(db_name, timeout=60.0)
            c = conn.cursor()
            c.execute('''INSERT OR REPLACE INTO polyconv_states (n, d, b, k_start, logZ, prev_logZn, small_count, cumulative_capacity)
                         VALUES (?, ?, ?, ?, ?, ?, ?, ?)''', 
                      (n, d, b, k_start, logZ.tobytes(), prev_logZn, small_count, cumulative_capacity))
            conn.commit()
            conn.close()
            break
        except sqlite3.OperationalError:
            time.sleep(2)

def load_polyconv_state(n, d, b, db_name='recursions_polyconv.db'):
    for _ in range(5):
        try:
            conn = sqlite3.connect(db_name, timeout=60.0)
            c = conn.cursor()
            c.execute('''SELECT k_start, logZ, prev_logZn, small_count, cumulative_capacity FROM polyconv_states WHERE n=? AND b=? AND d=?''', (n, b, d))
            row = c.fetchone()
            conn.close()
            if row is None:
                return None
            k_start = row[0]
            logZ = np.frombuffer(row[1], dtype=np.float64).copy()
            prev_logZn = row[2]
            small_count = row[3]
            cumulative_capacity = row[4]
            return k_start, logZ, prev_logZn, small_count, cumulative_capacity
        except sqlite3.OperationalError:
            time.sleep(2)
    return None

@numba.njit
def logaddexp(a, b):
    if a == -np.inf:
        return b
    if b == -np.inf:
        return a
    if a > b:
        return a + math.log1p(math.exp(b - a))
    else:
        return b + math.log1p(math.exp(a - b))

@numba.njit
def level_degeneracy_numba(d, k):
    if k == 0: return 1.0
    if d == 1: return 1.0
    res = 1.0
    for i in range(1, k + 1):
        res = res * (d + k - i) / i
    return res

@numba.njit
def log_binom(n, k):
    if k < 0 or k > n:
        return -np.inf
    return math.lgamma(n + 1) - math.lgamma(k + 1) - math.lgamma(n - k + 1)

@numba.njit
def _polyconv_fixed_chunk(k_start, k_end, logZ, d, logb, n):
    for k in range(k_start, k_end):
        g_float = level_degeneracy_numba(d, k)
        mmax = int(min(g_float, float(n)))
        logwk = k * logb
        
        logC = np.zeros(mmax + 1, dtype=np.float64)
        for m in range(1, mmax + 1):
            logC[m] = log_binom(g_float, m) + m * logwk
            
        old = logZ.copy()
        for r in range(n + 1):
            s = -np.inf
            upper = min(r, mmax)
            for m in range(upper + 1):
                if old[r - m] != -np.inf:
                    s = logaddexp(s, logC[m] + old[r - m])
            logZ[r] = s
    return logZ, k_end

@numba.njit
def _polyconv_adaptive_chunk(k_start, k_end, logZ, d, logb, n, tol, consecutive_small, cumulative_capacity, prev_logZn, small_count):
    stop_reached = False
    for k in range(k_start, k_end):
        g_float = level_degeneracy_numba(d, k)
        mmax = int(min(g_float, float(n)))
        logwk = k * logb
        
        logC = np.zeros(mmax + 1, dtype=np.float64)
        for m in range(1, mmax + 1):
            logC[m] = log_binom(g_float, m) + m * logwk
            
        old = logZ.copy()
        for r in range(n + 1):
            s = -np.inf
            upper = min(r, mmax)
            for m in range(upper + 1):
                if old[r - m] != -np.inf:
                    s = logaddexp(s, logC[m] + old[r - m])
            logZ[r] = s

        cumulative_capacity += g_float

        if cumulative_capacity >= n:
            curr = logZ[n]
            if prev_logZn != -np.inf and curr != -np.inf:
                delta_log = abs(curr - prev_logZn)
                if delta_log <= tol:
                    small_count += 1
                else:
                    small_count = 0
            prev_logZn = curr

            if small_count >= consecutive_small:
                stop_reached = True
                k_start = k + 1
                break
    else:
        k_start = k_end
        
    return logZ, cumulative_capacity, prev_logZn, small_count, stop_reached, k_start

def fermion_logZ_numeric(tau, N, n, d, max_shell=None, tol=1e-4,
                         consecutive_small=8, safety_cap=100000,
                         return_all=False):
    if N <= 0:
        raise ValueError("N must be a positive integer.")
    if n < 0:
        raise ValueError("n must be a nonnegative integer.")
    if d <= 0:
        raise ValueError("d must be a positive integer.")
    if tau < 0:
        raise ValueError("tau must be nonnegative.")
    if max_shell is not None and max_shell < 0:
        raise ValueError("max_shell must be None or a nonnegative integer.")

    epsilon = tau / N
    zeta = 1.0 + 0.5 * epsilon * epsilon
    u = math.log(zeta + math.sqrt(zeta * zeta - 1.0))
    b = math.exp(-N * u)
    logb = -N * u

    if n == 0:
        logZ = [-math.inf] * (n + 1)
        logZ[0] = 0.0
        info = {
            "epsilon": epsilon,
            "zeta": zeta,
            "u": u,
            "b": b,
            "logb": logb,
            "shells_used": 0,
        }
        return (0.0, logZ, info) if return_all else 0.0

    init_polyconv_db()
    state = load_polyconv_state(n, d, b)
    if state is not None:
        k_start, logZ, prev_logZn, small_count, cumulative_capacity = state
    else:
        k_start = 0
        logZ = np.full(n + 1, -np.inf, dtype=np.float64)
        logZ[0] = 0.0
        prev_logZn = -np.inf
        small_count = 0
        cumulative_capacity = 0.0

    chunk_size = 50000
    stop_reached = False

    if max_shell is not None:
        target_cap = max_shell + 1
        is_adaptive = False
    else:
        target_cap = safety_cap + 1
        is_adaptive = True

    while k_start < target_cap and not stop_reached:
        k_end = min(k_start + chunk_size, target_cap)
        if is_adaptive:
            logZ, cumulative_capacity, prev_logZn, small_count, stop_reached, k_start = _polyconv_adaptive_chunk(
                k_start, k_end, logZ, d, logb, n, tol, consecutive_small, cumulative_capacity, prev_logZn, small_count
            )
        else:
            logZ, k_start = _polyconv_fixed_chunk(k_start, k_end, logZ, d, logb, n)
            
        save_polyconv_state(n, d, b, k_start, logZ, prev_logZn, small_count, cumulative_capacity)

    if is_adaptive and not stop_reached:
        raise RuntimeError(
            "Adaptive truncation did not converge. "
            "Try increasing safety_cap or set max_shell explicitly."
        )

    shells_used = k_start

    logZ_list = [float(x) if x != -np.inf else -math.inf for x in logZ]

    if return_all:
        info = {
            "epsilon": epsilon,
            "zeta": zeta,
            "u": u,
            "b": b,
            "logb": logb,
            "shells_used": shells_used,
        }
        return logZ_list[n], logZ_list, info

    return logZ_list[n]


import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from joblib import Parallel, delayed

def make_tau_grid(tau_start, tau_end, tau_step):
    if tau_step <= 0:
        raise ValueError("tau_step must be positive.")
    if tau_end < tau_start:
        raise ValueError("tau_end must be >= tau_start.")

    taus = []
    tau = tau_start
    while tau <= tau_end + 1e-15:
        taus.append(round(tau, 12))
        tau += tau_step
    return taus

def finite_difference_user_sign(x, y):
    if len(x) != len(y):
        raise ValueError("x and y must have the same length.")
    if len(x) < 2:
        raise ValueError("Need at least 2 points.")

    x_fd = x[:-1]
    y_fd = []
    for i in range(len(x) - 1):
        dx = x[i + 1] - x[i]
        y_fd.append((y[i] - y[i + 1]) / dx)
    return np.array(x_fd), np.array(y_fd)

def compute_logZ_worker(tau, N, n, d, logZ_kwargs):
    return fermion_logZ_numeric(tau=tau, N=N, n=n, d=d, **logZ_kwargs)

def compute_logZ_worker_all(tau, N, n, d, logZ_kwargs):
    logZ_kwargs['return_all'] = True
    _, logZ_list, _ = fermion_logZ_numeric(tau=tau, N=N, n=n, d=d, **logZ_kwargs)
    return logZ_list

def plot_logZ_and_fd_multiN(tau_start, tau_end, tau_step, n, d, N_list,
                            show_logZ=True, show_fd=True, store_all_n=False, save_filename="plot_data", **logZ_kwargs):
    taus = np.array(make_tau_grid(tau_start, tau_end, tau_step))
    results = {}

    if store_all_n:
        results_fd_all_n = []
        results_logZ_all_n = []

    for N in N_list:
        if store_all_n:
            logZ_kwargs['return_all'] = True
            logZ_vals = Parallel(n_jobs=-1)(
                delayed(compute_logZ_worker_all)(tau, N, n, d, logZ_kwargs)
                for tau in taus
            )
            # logZ_vals is list of length num_taus, each element is list of length n+1.
            logZ_vals_arr = np.array(logZ_vals).T # shape: (n+1, num_taus)
            
            fd_vals_arr = []
            for i in range(n + 1):
                tau_fd, fd_vals = finite_difference_user_sign(taus, logZ_vals_arr[i])
                fd_vals_arr.append(fd_vals)
            fd_vals_arr = np.array(fd_vals_arr) # shape: (n+1, num_taus - 1)
            
            results_logZ_all_n.append(logZ_vals_arr)
            results_fd_all_n.append(fd_vals_arr)
            
            results[N] = {
                "taus": taus,
                "logZ": logZ_vals_arr[n],
                "tau_fd": tau_fd,
                "fd": fd_vals_arr[n],
            }
        else:
            logZ_kwargs['return_all'] = False
            logZ_vals = Parallel(n_jobs=-1)(
                delayed(compute_logZ_worker)(tau, N, n, d, logZ_kwargs)
                for tau in taus
            )
            tau_fd, fd_vals = finite_difference_user_sign(taus, logZ_vals)
    
            results[N] = {
                "taus": taus,
                "logZ": logZ_vals,
                "tau_fd": tau_fd,
                "fd": fd_vals,
            }

    if show_logZ:
        plt.figure(figsize=(8, 5))
        for N in N_list:
            plt.plot(results[N]["taus"], results[N]["logZ"], marker=".", label=f"N={N}")
        plt.xlabel("tau")
        plt.ylabel("log Z_n")
        plt.title(f"log Z_n vs tau (n={n}, d={d})")
        plt.grid(True)
        plt.legend()
        plt.show()

    if show_fd:
        plt.figure(figsize=(8, 5))
        for N in N_list:
            plt.plot(results[N]["tau_fd"], results[N]["fd"], marker=".", label=f"N={N}")
        plt.xlabel("tau")
        plt.ylabel("(log Z(tau) - log Z(tau + dt)) / dt")
        plt.title(f"Finite difference vs tau (n={n}, d={d})")
        plt.grid(True)
        plt.legend()
        plt.show()

    if store_all_n:
        arr_logZ = np.array(results_logZ_all_n)[:, 1:, :].transpose(1, 0, 2)
        arr_fd = np.array(results_fd_all_n)[:, 1:, :].transpose(1, 0, 2)
        
        np.save(f"{save_filename}_logZ_all_n.npy", arr_logZ)
        np.save(f"{save_filename}_fd_all_n.npy", arr_fd)
        np.save(f"{save_filename}_taus.npy", taus)
        np.save(f"{save_filename}_tau_fd.npy", tau_fd)
        np.save(f"{save_filename}_beads.npy", np.array(N_list))
        
        print(f"Data successfully saved to '{save_filename}_fd_all_n.npy' (shape {arr_fd.shape})")
        print(f"Axes mapping: n (axis 0), beads (axis 1) -> '{save_filename}_beads.npy', taus -> '{save_filename}_tau_fd.npy'")

    return results

res = plot_logZ_and_fd_multiN(tau_start=0.1, tau_end=1.0, tau_step=0.5, n=7, d=2, N_list=[1, 4], show_logZ=False, show_fd=False, store_all_n=True)
