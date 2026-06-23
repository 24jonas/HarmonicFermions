"""
Runtime scaling benchmarks for the PolyConv method.

Provides functions for measuring and fitting runtime as a function of:
  - n (particle number)
  - tau (imaginary time)
  - N (bead count)

Used by PolyConv-ArbitraryPotential and PolyConv-SphericalPotential.
"""

import time

import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit

# ============================================================
# Scaling vs n
# ============================================================


def scaling_model(n, a, b, c):
    """Power-law model: t(n) = a * n^b + c."""
    return a * (n**b) + c


def fit_runtime_scaling(
    tau,
    N,
    d,
    logZ_func,
    n_start=0,
    n_end=2000,
    n_step=25,
    repeats=3,
    use_median=True,
    **logZ_kwargs,
):
    """Measure and fit runtime vs n."""
    n_vals = np.arange(n_start, n_end + 1, n_step)
    runtimes = []

    for n_val in n_vals:
        samples = []
        for _ in range(repeats):
            t0 = time.perf_counter()
            logZ_func(tau=tau, N=N, n=int(n_val), d=d, **logZ_kwargs)
            t1 = time.perf_counter()
            samples.append(t1 - t0)
        runtime = float(np.median(samples) if use_median else np.mean(samples))
        runtimes.append(runtime)

    runtimes = np.array(runtimes, dtype=float)
    c0 = float(runtimes.min())
    a0 = max(float(runtimes.max() - c0), 1e-12)
    b0 = 1.0
    popt, pcov = curve_fit(
        scaling_model, n_vals.astype(float), runtimes, p0=(a0, b0, c0), maxfev=20000
    )
    return {
        "n_vals": n_vals,
        "runtimes": runtimes,
        "params": {"a": popt[0], "b": popt[1], "c": popt[2]},
        "cov": pcov,
    }


def plot_runtime_scaling_fit(
    tau,
    N,
    d,
    logZ_func,
    n_start=0,
    n_end=2000,
    n_step=25,
    repeats=3,
    use_median=True,
    **logZ_kwargs,
):
    """Measure, fit, and plot runtime vs n."""
    result = fit_runtime_scaling(
        tau=tau,
        N=N,
        d=d,
        logZ_func=logZ_func,
        n_start=n_start,
        n_end=n_end,
        n_step=n_step,
        repeats=repeats,
        use_median=use_median,
        **logZ_kwargs,
    )
    n_vals = result["n_vals"]
    runtimes = result["runtimes"]
    a, b, c = result["params"]["a"], result["params"]["b"], result["params"]["c"]

    n_fit = np.linspace(n_vals.min(), n_vals.max(), 400)
    t_fit = scaling_model(n_fit, a, b, c)

    plt.figure(figsize=(8, 5))
    plt.plot(n_vals, runtimes, "o", label="measured runtime")
    plt.plot(n_fit, t_fit, label=f"fit: a*n^b + c\n a={a:.6g}, b={b:.6g}, c={c:.6g}")
    plt.xlabel("n")
    plt.ylabel("runtime (seconds)")
    plt.title(f"Runtime scaling at tau={tau}, N={N}, d={d}")
    plt.grid(True)
    plt.legend()
    plt.show()
    print(f"a = {a}\nb = {b}\nc = {c}")
    return result


# ============================================================
# Scaling vs tau
# ============================================================


def scaling_model_tau_negpow(tau, a, b, c):
    """Negative power-law model: t(tau) = a * tau^(-b) + c."""
    return a * (tau ** (-b)) + c


def highest_shell_zero_temp(n, d=2):
    """Return the highest occupied shell at T=0 for n fermions in d dimensions."""
    total = 0
    k = 0
    while total < n:
        if d == 1:
            g_k = 1
        else:
            g_k = 1
            for i in range(1, k + 1):
                g_k = g_k * (d + k - i) // i
        total += g_k
        k += 1
    return k


def fit_runtime_vs_tau_negpow(
    tau_start,
    tau_end,
    tau_step,
    N,
    n,
    d,
    logZ_func,
    repeats=3,
    use_median=True,
    **logZ_kwargs,
):
    """Measure and fit runtime vs tau."""
    tau_vals = np.arange(tau_start, tau_end + 0.5 * tau_step, tau_step, dtype=float)
    runtimes = []
    shell_counts = []

    for tau in tau_vals:
        samples = []
        shell_used = None
        for _ in range(repeats):
            t0 = time.perf_counter()
            _, _, info = logZ_func(
                tau=float(tau), N=N, n=n, d=d, return_all=True, **logZ_kwargs
            )
            t1 = time.perf_counter()
            samples.append(t1 - t0)
            if shell_used is None:
                shell_used = info["shells_used"]
        runtime = float(np.median(samples) if use_median else np.mean(samples))
        runtimes.append(runtime)
        shell_counts.append(shell_used)

    runtimes = np.array(runtimes, dtype=float)
    shell_counts = np.array(shell_counts, dtype=int)
    c0 = float(runtimes.min())
    a0 = max(float(runtimes.max() - c0), 1e-12)
    popt, pcov = curve_fit(
        scaling_model_tau_negpow, tau_vals, runtimes, p0=(a0, 1.0, c0), maxfev=20000
    )
    return {
        "tau_vals": tau_vals,
        "runtimes": runtimes,
        "shell_counts": shell_counts,
        "params": {"a": popt[0], "b": popt[1], "c": popt[2]},
        "cov": pcov,
    }


def plot_runtime_vs_tau_negpow_fit(
    tau_start,
    tau_end,
    tau_step,
    N,
    n,
    d,
    logZ_func,
    repeats=3,
    use_median=True,
    print_shells=True,
    **logZ_kwargs,
):
    """Measure, fit, and plot runtime vs tau."""
    result = fit_runtime_vs_tau_negpow(
        tau_start=tau_start,
        tau_end=tau_end,
        tau_step=tau_step,
        N=N,
        n=n,
        d=d,
        logZ_func=logZ_func,
        repeats=repeats,
        use_median=use_median,
        **logZ_kwargs,
    )
    tau_vals = result["tau_vals"]
    runtimes = result["runtimes"]
    shell_counts = result["shell_counts"]
    a, b, c = result["params"]["a"], result["params"]["b"], result["params"]["c"]

    tau_fit = np.linspace(tau_vals.min(), tau_vals.max(), 400)
    t_fit = scaling_model_tau_negpow(tau_fit, a, b, c)

    plt.figure(figsize=(8, 5))
    plt.plot(tau_vals, runtimes, "o", label="measured runtime")
    plt.plot(
        tau_fit, t_fit, label=f"fit: a*tau^(-b)+c\n a={a:.6g}, b={b:.6g}, c={c:.6g}"
    )
    plt.xlabel("tau")
    plt.ylabel("runtime (seconds)")
    plt.title(f"Runtime vs tau (n={n}, N={N}, d={d})")
    plt.grid(True)
    plt.legend()
    plt.show()

    if print_shells:
        print("Shells used per tau:")
        for t, s in zip(tau_vals, shell_counts):
            print(f"  tau={t:.4f} -> shells={s}")
    print(f"a = {a}\nb = {b}\nc = {c}")
    return result


# ============================================================
# Scaling vs N (bead count)
# ============================================================


def scaling_model_N_negpow(N, a, b, c):
    """Negative power-law model: t(N) = a * N^(-b) + c."""
    return a * (N ** (-b)) + c


def fit_runtime_vs_N_negpow(
    N_start,
    N_end,
    N_step,
    tau,
    n,
    d,
    logZ_func,
    repeats=3,
    use_median=True,
    **logZ_kwargs,
):
    """Measure and fit runtime vs N."""
    N_vals = np.arange(N_start, N_end + 0.5 * N_step, N_step, dtype=int)
    runtimes = []
    shell_counts = []

    for N in N_vals:
        samples = []
        shell_used = None
        for _ in range(repeats):
            t0 = time.perf_counter()
            _, _, info = logZ_func(
                tau=tau, N=int(N), n=n, d=d, return_all=True, **logZ_kwargs
            )
            t1 = time.perf_counter()
            samples.append(t1 - t0)
            if shell_used is None:
                shell_used = info["shells_used"]
        runtime = float(np.median(samples) if use_median else np.mean(samples))
        runtimes.append(runtime)
        shell_counts.append(shell_used)

    runtimes = np.array(runtimes, dtype=float)
    shell_counts = np.array(shell_counts, dtype=int)
    c0 = float(runtimes.min())
    a0 = max(float(runtimes.max() - c0), 1e-12)
    popt, pcov = curve_fit(
        scaling_model_N_negpow,
        N_vals.astype(float),
        runtimes,
        p0=(a0, 1.0, c0),
        maxfev=20000,
    )
    return {
        "N_vals": N_vals,
        "runtimes": runtimes,
        "shell_counts": shell_counts,
        "params": {"a": popt[0], "b": popt[1], "c": popt[2]},
        "cov": pcov,
    }


def plot_runtime_vs_N_negpow_fit(
    N_start,
    N_end,
    N_step,
    tau,
    n,
    d,
    logZ_func,
    repeats=3,
    use_median=True,
    print_shells=True,
    **logZ_kwargs,
):
    """Measure, fit, and plot runtime vs N."""
    result = fit_runtime_vs_N_negpow(
        N_start=N_start,
        N_end=N_end,
        N_step=N_step,
        tau=tau,
        n=n,
        d=d,
        logZ_func=logZ_func,
        repeats=repeats,
        use_median=use_median,
        **logZ_kwargs,
    )
    N_vals = result["N_vals"]
    runtimes = result["runtimes"]
    shell_counts = result["shell_counts"]
    a, b, c = result["params"]["a"], result["params"]["b"], result["params"]["c"]

    N_fit = np.linspace(N_vals.min(), N_vals.max(), 400)
    t_fit = scaling_model_N_negpow(N_fit, a, b, c)

    plt.figure(figsize=(8, 5))
    plt.plot(N_vals, runtimes, "o", label="measured runtime")
    plt.plot(N_fit, t_fit, label=f"fit: a*N^(-b)+c\n a={a:.6g}, b={b:.6g}, c={c:.6g}")
    plt.xlabel("N (beads)")
    plt.ylabel("runtime (seconds)")
    plt.title(f"Runtime vs N (n={n}, tau={tau}, d={d})")
    plt.grid(True)
    plt.legend()
    plt.show()

    if print_shells:
        print("Shells used per N:")
        for nv, s in zip(N_vals, shell_counts):
            print(f"  N={nv} -> shells={s}")
    print(f"a = {a}\nb = {b}\nc = {c}")
    return result
