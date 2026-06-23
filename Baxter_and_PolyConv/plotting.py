"""
Shared plotting and analysis utilities.

Contains helper functions for constructing tau grids, computing finite
differences, heat capacities, and producing logZ / energy plots.
"""

import math
import numpy as np


def make_tau_grid(tau_start, tau_end, tau_step):
    """Create an evenly spaced grid of tau values."""
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
    """
    Compute (y_i - y_{i+1}) / (x_{i+1} - x_i).
    Appropriate for E = -d(logZ)/d(tau).
    """
    if len(x) != len(y):
        raise ValueError("x and y must have the same length.")
    if len(x) < 2:
        raise ValueError("Need at least 2 points.")
    x_fd = x[:-1]
    y_fd = []
    for i in range(len(x) - 1):
        dx = x[i + 1] - x[i]
        y_fd.append((y[i] - y[i + 1]) / dx)
    return x_fd, y_fd


def compute_heat_capacity(tau_vals, E_vals):
    """Compute C_v = -tau^2 * dE/dtau from energy vs tau data."""
    tau_vals = np.asarray(tau_vals)
    E_vals = np.asarray(E_vals)
    dtau = np.diff(tau_vals)
    dE = np.diff(E_vals)
    tau_mid = (tau_vals[:-1] + tau_vals[1:]) / 2.0
    Cv = -(tau_mid**2) * (dE / dtau)
    return Cv


def plot_logZ_and_fd_multiN(tau_start, tau_end, tau_step, n, d, N_list,
                            show_logZ=True, show_fd=True,
                            logZ_func=None, **logZ_kwargs):
    """Plot log Z_n and/or its finite difference across multiple bead counts N."""
    import matplotlib.pyplot as plt

    if logZ_func is None:
        raise ValueError("logZ_func must be provided.")

    taus = make_tau_grid(tau_start, tau_end, tau_step)
    results = {}

    for N in N_list:
        logZ_vals = [logZ_func(tau=tau, N=N, n=n, d=d, **logZ_kwargs) for tau in taus]
        tau_fd, fd_vals = finite_difference_user_sign(taus, logZ_vals)
        results[N] = {"taus": taus, "logZ": logZ_vals, "tau_fd": tau_fd, "fd": fd_vals}

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

    return results
