"""
Physical quantities for non-interacting fermions in harmonic traps.

Provides exact ground-state energies (2D and d-dimensional), the
Thomas-Fermi (semiclassical) approximation for thermodynamic quantities,
and shell-filling utilities.
"""



import numpy as np
from scipy.integrate import quad
from scipy.optimize import root_scalar

# ============================================================
# Exact Ground-State Energies
# ============================================================


def energy_2d(n, w, verbose=False):
    """
    Exact energy for n non-interacting fermions in a 2D harmonic trap
    with relative frequency w.

    The energy is decomposed as:
        E(n, w) = E_CM + w * E_rel
    where E_CM = 1.0 (the single-particle center-of-mass energy in 2D).

    Parameters
    ----------
    n : int
        Number of fermions.
    w : float
        Effective relative frequency.
    verbose : bool
        If True, print the highest occupied shell.

    Returns
    -------
    float
        Total energy.
    """
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

    # Scale the relative energy component by w
    # For a 2D harmonic oscillator, the 1-particle CM energy is exactly 1.0
    e_1 = 1.0
    true_energy_w = e_1 + w * (base_energy - e_1)

    return true_energy_w


def energy_nd(n, d, w, verbose=False):
    """
    Exact energy for n non-interacting fermions in a d-dimensional
    harmonic trap with frequency w.

    Shell-filling with degeneracy g_k = C(d + k - 1, k).

    Parameters
    ----------
    n : int
        Number of fermions.
    d : int
        Spatial dimension.
    w : float
        Harmonic frequency.
    verbose : bool
        If True, print the highest occupied shell.

    Returns
    -------
    float
        Total energy E = w * (sum of occupied level energies + n * d/2).
    """
    total_particles = 0
    total_energy = 0
    k = 0
    while total_particles < n:
        if d == 1:
            g_k = 1
        else:
            g_k = 1
            for i in range(1, k + 1):
                g_k = g_k * (d + k - i) // i

        take = min(g_k, n - total_particles)
        total_particles += take
        total_energy += take * k
        k += 1

    if verbose:
        print("Highest shell:", k)

    return w * (total_energy + n * d / 2.0)


def highest_occupied_shell(n, d=2):
    """
    Return the index of the highest occupied energy shell when filling
    n fermions into a d-dimensional harmonic oscillator at zero temperature.

    Parameters
    ----------
    n : int
        Number of fermions.
    d : int
        Spatial dimension (default 2).

    Returns
    -------
    int
        Highest shell index (0-based).
    """
    total_particles = 0
    k = 0
    while total_particles < n:
        if d == 1:
            g_k = 1
        else:
            g_k = 1
            for i in range(1, k + 1):
                g_k = g_k * (d + k - i) // i
        total_particles += g_k
        k += 1
    return k


# ============================================================
# Thomas-Fermi (Semiclassical) Approximation
# ============================================================


def N_int(mu, T, w=1.0):
    """
    Thomas-Fermi particle number: integral of the Fermi-Dirac density of
    states in a 2D harmonic trap.

    Parameters
    ----------
    mu : float
        Chemical potential.
    T : float
        Temperature.
    w : float
        Harmonic frequency.

    Returns
    -------
    float
        Particle number N(mu, T).
    """

    def integrand(e):
        x = (e - mu) / T
        if x > 100:
            return 0.0
        return (e / w**2) / (np.exp(x) + 1.0)

    return quad(integrand, 0, max(0, mu) + 40 * T)[0]


def E_int(mu, T, w=1.0):
    """
    Thomas-Fermi total energy: integral of the Fermi-Dirac density of
    states weighted by energy in a 2D harmonic trap.

    Parameters
    ----------
    mu : float
        Chemical potential.
    T : float
        Temperature.
    w : float
        Harmonic frequency.

    Returns
    -------
    float
        Total energy E(mu, T).
    """

    def integrand(e):
        x = (e - mu) / T
        if x > 100:
            return 0.0
        return (e**2 / w**2) / (np.exp(x) + 1.0)

    return quad(integrand, 0, max(0, mu) + 40 * T)[0]


def get_mu_TF(n_target, T_fixed, w=1.0):
    """
    Solve for the Thomas-Fermi chemical potential at fixed T and n.

    Parameters
    ----------
    n_target : int or float
        Target particle number.
    T_fixed : float
        Temperature.
    w : float
        Harmonic frequency.

    Returns
    -------
    float
        Chemical potential mu.
    """
    mu0 = np.sqrt(2 * n_target) * w

    def obj(mu):
        return N_int(mu, T_fixed, w) - n_target

    res = root_scalar(obj, bracket=[-20 * T_fixed, mu0 + 20 * T_fixed], method="brentq")
    return res.root


def get_E_TF(n_target, T_fixed, w=1.0):
    """
    Compute the Thomas-Fermi total energy at fixed T and n.

    Parameters
    ----------
    n_target : int or float
        Target particle number.
    T_fixed : float
        Temperature.
    w : float
        Harmonic frequency.

    Returns
    -------
    float
        Total energy.
    """
    mu = get_mu_TF(n_target, T_fixed, w)
    return E_int(mu, T_fixed, w)


def get_Cv_TF(n_target, T_fixed, w=1.0, dT=0.01):
    """
    Compute the Thomas-Fermi heat capacity at fixed T and n via finite
    difference.

    Parameters
    ----------
    n_target : int or float
        Target particle number.
    T_fixed : float
        Temperature.
    w : float
        Harmonic frequency.
    dT : float
        Temperature step for finite difference.

    Returns
    -------
    float
        Heat capacity C_v.
    """
    T2 = T_fixed + dT
    E1 = get_E_TF(n_target, T_fixed, w)
    E2 = get_E_TF(n_target, T2, w)
    return (E2 - E1) / dT
