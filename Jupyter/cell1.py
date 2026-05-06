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
            logZ_vals_arr = np.array(logZ_vals).T 
            
            fd_vals_arr = []
            for i in range(n + 1):
                tau_fd, fd_vals = finite_difference_user_sign(taus, logZ_vals_arr[i])
                fd_vals_arr.append(fd_vals)
            fd_vals_arr = np.array(fd_vals_arr)
            
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
    else:
        df = pd.DataFrame(results[N_list[0]])
        df.to_csv(f"{save_filename}_results.csv", index=False)

    return results
