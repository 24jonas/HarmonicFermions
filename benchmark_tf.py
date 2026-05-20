import numpy as np
import time
from scipy.integrate import quad
from scipy.optimize import root_scalar

w_val = 0.5
def N_int(mu, T):
    def integrand(e):
        x = (e - mu)/T
        if x > 100: return 0.0
        return (e / w_val**2) / (np.exp(x) + 1.0)
    return quad(integrand, 0, max(0, mu) + 40*T)[0]

def get_mu_TF(n_target, T_fixed):
    mu0 = np.sqrt(2 * n_target) * w_val
    def obj1(mu): return N_int(mu, T_fixed) - n_target
    res1 = root_scalar(obj1, bracket=[-20*T_fixed, mu0 + 20*T_fixed], method='brentq')
    return res1.root

start = time.time()
for n in range(1, 1001):
    get_mu_TF(n, 0.5)
end = time.time()
print(f"Time for 1000 points: {end - start:.2f} seconds")
