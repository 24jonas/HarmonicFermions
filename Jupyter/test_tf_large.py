import numpy as np
import scipy.integrate as integrate
import scipy.optimize as optimize
from scipy.special import expit

def F1(eta):
    if eta > 0:
        return integrate.quad(lambda x: x * expit(eta - x), 0, eta, limit=200)[0] + \
               integrate.quad(lambda x: x * expit(eta - x), eta, eta+50, limit=200)[0]
    else:
        return integrate.quad(lambda x: x * expit(eta - x), 0, 50, limit=200)[0]

def F2(eta):
    if eta > 0:
        return integrate.quad(lambda x: 0.5 * x**2 * expit(eta - x), 0, eta, limit=200)[0] + \
               integrate.quad(lambda x: 0.5 * x**2 * expit(eta - x), eta, eta+50, limit=200)[0]
    else:
        return integrate.quad(lambda x: 0.5 * x**2 * expit(eta - x), 0, 50, limit=200)[0]

def solve_eta(n_val, T):
    target = n_val / (T**2)
    eta_guess = np.sqrt(2 * target)
    res = optimize.root_scalar(lambda eta: F1(eta) - target, x0=eta_guess, x1=eta_guess*1.01)
    return res.root

T = 0.066
n_val = 10000
eta = solve_eta(n_val, T)
print(f"eta={eta}")

