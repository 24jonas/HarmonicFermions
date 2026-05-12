import numpy as np
from scipy.integrate import quad
from scipy.optimize import root_scalar

def N_int(mu, T):
    def integrand(e):
        x = (e - mu)/T
        if x > 100: return 0.0
        return e / (np.exp(x) + 1.0)
    return quad(integrand, 0, max(0, mu) + 40*T)[0]

N = 10000
T = 10.0
mu0 = np.sqrt(2*N)
def obj(mu):
    return N_int(mu, T) - N
res = root_scalar(obj, bracket=[-20*T, mu0 + 20*T], method='brentq')
print(res)
