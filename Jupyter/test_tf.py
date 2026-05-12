import numpy as np
from scipy.integrate import quad
from scipy.optimize import root_scalar

def N_int(mu, T):
    # n = -T^2 * Li2(-exp(mu/T))
    # But for numerical stability at large mu/T, integration is safe
    def integrand(e):
        # to avoid overflow:
        x = (e - mu)/T
        if x > 100: return 0.0
        return e / (np.exp(x) + 1.0)
    return quad(integrand, 0, mu + 40*T)[0]

def E_int(mu, T):
    def integrand(e):
        x = (e - mu)/T
        if x > 100: return 0.0
        return (e**2) / (np.exp(x) + 1.0)
    return quad(integrand, 0, mu + 40*T)[0]

def calc_TF(N, T_array):
    E_TF = np.zeros_like(T_array)
    mu_TF = np.zeros_like(T_array)
    for i, T in enumerate(T_array):
        # find mu
        # at T=0, mu = sqrt(2N)
        mu0 = np.sqrt(2*N)
        def obj(mu):
            return N_int(mu, T) - N
        res = root_scalar(obj, bracket=[0, mu0 + 10*T], method='brentq')
        mu_TF[i] = res.root
        E_TF[i] = E_int(mu_TF[i], T)
        
    return E_TF, mu_TF

T_array = np.linspace(0.1, 10, 20)
E_TF, mu_TF = calc_TF(10000, T_array)
print("T:", T_array[:5])
print("E_TF/N:", (E_TF/10000)[:5])
