import numpy as np
import scipy.linalg as la
import scipy.special

epsilon = 2.5 / 4.0
r_min = 1e-6
r_max = 5.0
num_grid = 800

r = np.linspace(r_min, r_max, num_grid)
dr = (r_max - r_min) / (num_grid - 1)

R, Rp = np.meshgrid(r, r)

def V(x):
    return 0.5 * x**2

z = R * Rp / epsilon
l = 0

rho_l = (R * Rp / epsilon) * np.exp(-((R - Rp)**2) / (2 * epsilon)) * scipy.special.ive(l + 0.5, z)
T_l = np.exp(-epsilon * V(R) / 2) * rho_l * np.exp(-epsilon * V(Rp) / 2) * dr

eigvals = la.eigvalsh(T_l)
print(np.sort(eigvals)[::-1][:5])
