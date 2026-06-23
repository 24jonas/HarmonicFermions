"""
Polynomial Convolution (PolyConv) method — General partition function
computation for fermions with arbitrary single-particle spectra.

Computes log Z_n by convolving per-shell generating functions in
log-space. Accepts custom degeneracy (gk_fn) and log-weight (logwk_fn)
functions, making it applicable to arbitrary potentials (not just
harmonic).

For the harmonic-specific, numba-accelerated version with DB caching,
see polyconv_opt.py.
"""

import math


# ============================================================
# Log-space arithmetic utilities
# ============================================================

def logaddexp(a, b):
    """Numerically stable log(exp(a) + exp(b))."""
    if a == -math.inf:
        return b
    if b == -math.inf:
        return a
    if a > b:
        return a + math.log1p(math.exp(b - a))
    else:
        return b + math.log1p(math.exp(a - b))


# ============================================================
# Degeneracy functions
# ============================================================

def degeneracy_box(d, k):
    """
    Degeneracy of the k-th energy level for a d-dimensional
    isotropic harmonic oscillator (box/Cartesian spectrum).

    This is C(d + k - 1, k) = (d+k-1)! / (k! (d-1)!).
    """
    if k == 0:
        return 1
    if d == 1:
        return 1
    result = 1
    for i in range(1, k + 1):
        result = result * (d + k - i) // i
    return result


def level_degeneracy(d, k):
    """Alias for degeneracy_box (Cartesian harmonic oscillator)."""
    return degeneracy_box(d, k)


def log_binom(n, k):
    """Compute log(C(n, k)) using the log-gamma function."""
    if k < 0 or k > n:
        return -math.inf
    return math.lgamma(n + 1) - math.lgamma(k + 1) - math.lgamma(n - k + 1)


# ============================================================
# Core convolution
# ============================================================

def convolve_shell_logspace(logZ, g, logw, n):
    """
    Multiply current polynomial by (1 + w*t)^g in log-space,
    truncated to degree n.

    This convolves the running partition function with the contribution
    from a single energy shell of degeneracy g and log-weight logw.

    Parameters
    ----------
    logZ : list of float
        Current log-coefficients, length n+1.
    g : int
        Shell degeneracy.
    logw : float
        Log of the Boltzmann weight for this shell.
    n : int
        Maximum particle number (truncation degree).
    """
    if g < 0 or int(g) != g:
        raise ValueError(f"g_k must be a nonnegative integer, got {g}.")
    g = int(g)
    mmax = min(g, n)

    logC = [0.0] * (mmax + 1)
    for m in range(1, mmax + 1):
        logC[m] = log_binom(g, m) + m * logw

    for r in range(n, -1, -1):
        s = -math.inf
        upper = min(r, mmax)
        for m in range(upper + 1):
            if logZ[r - m] != -math.inf:
                s = logaddexp(s, logC[m] + logZ[r - m])
        logZ[r] = s


# ============================================================
# Main API
# ============================================================

def fermion_logZ_numeric(
    tau, N, n, d=None,
    max_shell=None, tol=1e-4, consecutive_small=8,
    safety_cap=100000,
    return_all=False,
    gk_fn=None, Ek_fn=None, logwk_fn=None,
    k_start=0,
):
    """
    Compute the fermionic canonical partition function log Z_n via
    polynomial convolution.

    This general implementation accepts custom per-shell degeneracy and
    weight functions, enabling use with arbitrary single-particle spectra
    (e.g., from numerical diagonalization).

    Parameters
    ----------
    tau : float
        Imaginary time (inverse temperature times bead count).
    N : int
        Number of Trotter beads.
    n : int
        Number of fermions.
    d : int or None
        Spatial dimension (required if gk_fn is not provided).
    max_shell : int or None
        Fixed maximum shell index. If None, uses adaptive truncation.
    tol : float
        Convergence tolerance for adaptive truncation.
    consecutive_small : int
        Number of consecutive small changes before declaring convergence.
    safety_cap : int
        Maximum shell index to prevent runaway computation.
    return_all : bool
        If True, return (logZ_n, logZ_list, info_dict).
    gk_fn : callable or None
        Function gk_fn(k) -> int giving degeneracy of shell k.
        Defaults to degeneracy_box(d, k).
    Ek_fn : callable or None
        Function Ek_fn(k) -> float giving energy of shell k.
        Used to compute logwk_fn if logwk_fn is not provided.
        Defaults to Ek_fn(k) = k.
    logwk_fn : callable or None
        Function logwk_fn(k, logb) -> float giving log-weight of shell k.
        Defaults to Ek_fn(k) * logb.
    k_start : int
        Starting shell index (for custom spectra).

    Returns
    -------
    float or tuple
        If return_all is False: logZ_n.
        If return_all is True: (logZ_n, logZ_list, info_dict).
    """
    if N <= 0:
        raise ValueError("N must be a positive integer.")
    if n < 0:
        raise ValueError("n must be a nonnegative integer.")
    if k_start < 0:
        raise ValueError("k_start must be nonnegative.")

    if gk_fn is None and d is None:
        raise ValueError("Provide d for the default degeneracy, or provide gk_fn.")
    if d is not None and d <= 0:
        raise ValueError("d must be a positive integer.")

    epsilon = tau / N
    zeta = 1.0 + 0.5 * epsilon * epsilon
    u = math.log(zeta + math.sqrt(zeta * zeta - 1.0))
    b = math.exp(-N * u)
    logb = -N * u  # log(b)

    # Defaults
    if gk_fn is None:
        gk_fn = lambda k: degeneracy_box(d, k)

    if logwk_fn is None:
        if Ek_fn is None:
            Ek_fn = lambda k: k
        logwk_fn = lambda k, logb: Ek_fn(k) * logb

    # logZ[m] = log coefficient of t^m in running polynomial
    logZ = [-math.inf] * (n + 1)
    logZ[0] = 0.0

    if n == 0:
        info = {
            "epsilon": epsilon,
            "zeta": zeta,
            "u": u,
            "b": b,
            "logb": logb,
            "shells_used": 0,
            "k_start": k_start,
        }
        return (0.0, logZ, info) if return_all else 0.0

    # Adaptive truncation state
    prev_logZn = -math.inf
    small_count = 0
    cumulative_capacity = 0.0

    if max_shell is not None:
        target = max_shell + 1
    else:
        target = safety_cap + 1

    k = k_start
    while k < target:
        g = gk_fn(k)
        logw = logwk_fn(k, logb)
        convolve_shell_logspace(logZ, g, logw, n)

        cumulative_capacity += g
        k += 1

        # Adaptive convergence check
        if max_shell is None and cumulative_capacity >= n:
            curr = logZ[n]
            if prev_logZn != -math.inf and curr != -math.inf:
                delta = abs(curr - prev_logZn)
                if delta <= tol:
                    small_count += 1
                else:
                    small_count = 0
            prev_logZn = curr

            if small_count >= consecutive_small:
                break

    shells_used = k

    if max_shell is None and small_count < consecutive_small:
        raise RuntimeError(
            "Adaptive truncation did not converge. "
            "Try increasing safety_cap or set max_shell explicitly."
        )

    # Apply zero-point energy shift
    logZ_out = list(logZ)
    for r in range(n + 1):
        if logZ_out[r] != -math.inf:
            if d is not None:
                logZ_out[r] += r * (d / 2.0) * logb
            # When d is None (custom spectrum), no shift is applied

    if return_all:
        info = {
            "epsilon": epsilon,
            "zeta": zeta,
            "u": u,
            "b": b,
            "logb": logb,
            "shells_used": shells_used,
            "k_start": k_start,
        }
        return logZ_out[n], logZ_out, info

    return logZ_out[n]
