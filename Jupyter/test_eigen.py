import numpy as np
import scipy.linalg as la
import math

def get_1d_eigenvalues(tau, N, L=2.0, num_grid=1000):
    epsilon = tau / N
    x = np.linspace(-L/2, L/2, num_grid)
    dx = L / (num_grid - 1)
    
    X, Y = np.meshgrid(x, x)
    # T matrix: T(x, y) * dx
    T = (1.0 / np.sqrt(2 * np.pi * epsilon)) * np.exp(-((X - Y)**2) / (2 * epsilon)) * dx
    
    # Eigenvalues
    eigvals = la.eigvalsh(T)
    # Filter positive
    eigvals = eigvals[eigvals > 1e-15]
    return np.sort(eigvals)[::-1] # descending

def get_2d_log_eigenvalues(tau, N, L=2.0, num_grid=500, max_states=400):
    eig1 = get_1d_eigenvalues(tau, N, L, num_grid)
    # 2D eigenvalues are products of 1D eigenvalues
    # log(lambda_2D) = log(lambda_1x) + log(lambda_1y)
    log_eig1 = np.log(eig1)
    
    log_eig2 = []
    for i in range(len(log_eig1)):
        for j in range(len(log_eig1)):
            if log_eig1[i] + log_eig1[j] > -100: # filter very small ones
                log_eig2.append(log_eig1[i] + log_eig1[j])
    
    log_eig2 = np.array(log_eig2)
    log_eig2.sort()
    log_eig2 = log_eig2[::-1] # descending
    return log_eig2[:max_states]

def logaddexp(a, b):
    if a == -math.inf: return b
    if b == -math.inf: return a
    if a > b: return a + math.log1p(math.exp(b - a))
    else: return b + math.log1p(math.exp(a - b))

def get_logZ_fermions(log_eig, N, n):
    # logw_k = N * log_lambda_k
    logw = N * log_eig
    logZ = [-math.inf] * (n + 1)
    logZ[0] = 0.0
    
    for lw in logw:
        old = logZ[:]
        for r in range(1, n + 1):
            if old[r - 1] != -math.inf:
                logZ[r] = logaddexp(old[r], old[r - 1] + lw)
    return logZ[n]

N_beads = 4
n_part = 5
L = 2.0

for tau in [2.0, 2.5]:
    log_eig = get_2d_log_eigenvalues(tau, N_beads, L, num_grid=1000)
    logZn = get_logZ_fermions(log_eig, N_beads, n_part)
    print(f"tau={tau}, logZ={logZn}")

logZ_20 = get_logZ_fermions(get_2d_log_eigenvalues(2.0, N_beads, L, 1000), N_beads, n_part)
logZ_25 = get_logZ_fermions(get_2d_log_eigenvalues(2.5, N_beads, L, 1000), N_beads, n_part)
print(f"E_GS for N=4: {(logZ_20 - logZ_25) / 0.5}")

