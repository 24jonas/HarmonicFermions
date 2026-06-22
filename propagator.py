"""
Primitive Approximation (PA) propagator functions and factor calculations.

These functions compute the propagator-specific quantities (zeta, lambda,
gamma, k1) and the thermodynamic/Hamiltonian correction factors used to
convert finite-difference estimators into physical energy estimators.

Used by both the Baxter method and PolyConv method notebooks.
"""

import math


# ============================================================
# PA Propagator Functions
# ============================================================

def p_funcs_zeta_1(e):
    """First PA propagator function: zeta_1(epsilon) = 1 + epsilon^2 / 2."""
    return 1.0 + (e**2) / 2.0


def p_funcs_lambda(e):
    """PA propagator lambda function (= 1 for primitive approximation)."""
    return 1.0


def p_funcs_gamma(e):
    """PA propagator gamma function: sqrt(1 + epsilon^2 / 4)."""
    return math.sqrt(1.0 + (e**2) / 4.0)


def p_funcs_k1(e):
    """PA propagator k1 function (= epsilon for primitive approximation)."""
    return e


# ============================================================
# Factor Calculations
# ============================================================

def factor_calc_T(lambda_val, gamma_val, w, lambda_val_s, gamma_val_s):
    """
    Compute the thermodynamic estimator correction factors.

    Returns
    -------
    tuple of (fT_reg, fT_star)
        Regular and starred thermodynamic factors.
    """
    return (lambda_val / gamma_val), (w * lambda_val_s / gamma_val_s)


def factor_calc_H(lambda_val, gamma_val, w, lambda_val_s, gamma_val_s):
    """
    Compute the Hamiltonian estimator correction factors.

    Returns
    -------
    tuple of (fH_reg, fH_star)
        Regular and starred Hamiltonian factors.
    """
    return 0.5 * (gamma_val + 1.0 / gamma_val), (w / 2.0) * (gamma_val_s + 1.0 / gamma_val_s)


# ============================================================
# Boltzmann factor from (tau, N, w)
# ============================================================

def get_b_val(tau_val, N_val, w_val):
    """
    Compute the Boltzmann factor b from imaginary time tau, bead count N,
    and frequency w.

    When N_val == 0, uses the exact exponential: b = exp(-w * |tau|).
    Otherwise, uses the PA propagator: b = exp(-N * acosh(zeta_1)).

    Parameters
    ----------
    tau_val : float
        Imaginary time.
    N_val : int
        Number of beads. Use 0 for the exact (continuous) limit.
    w_val : float
        Harmonic frequency.

    Returns
    -------
    float
        The Boltzmann factor b in (0, 1).
    """
    if N_val == 0:
        return math.exp(-w_val * abs(tau_val))
    eps = w_val * (tau_val / N_val)
    z1 = p_funcs_zeta_1(eps)
    u = math.acosh(z1) if z1 >= 1.0 else 0.0
    return math.exp(-N_val * u)
