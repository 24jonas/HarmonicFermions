import json

with open('JonasPolyConv-Copy3.ipynb', 'r') as f:
    nb = json.load(f)

eigen_code = """import numpy as np
import scipy.linalg as la

def get_1d_log_eigenvalues(tau, N, L=2.0, num_grid=1000):
    epsilon = tau / N
    x = np.linspace(-L/2, L/2, num_grid)
    dx = L / (num_grid - 1)
    
    X, Y = np.meshgrid(x, x)
    # Primitive approximation integral operator
    T = (1.0 / np.sqrt(2 * np.pi * epsilon)) * np.exp(-((X - Y)**2) / (2 * epsilon)) * dx
    
    eigvals = la.eigvalsh(T)
    eigvals = eigvals[eigvals > 1e-15] # filter numerical noise
    return np.log(np.sort(eigvals)[::-1])

def get_nd_log_eigenvalues(tau, N, d, L=2.0, num_grid=500, max_states=1000):
    log_eig1 = get_1d_log_eigenvalues(tau, N, L, num_grid)
    
    if d == 1:
        return log_eig1[:max_states]
    elif d == 2:
        # Sum of logs for products of eigenvalues
        log_eig2 = []
        for i in range(len(log_eig1)):
            for j in range(len(log_eig1)):
                if log_eig1[i] + log_eig1[j] > -100:
                    log_eig2.append(log_eig1[i] + log_eig1[j])
        log_eig2 = np.array(log_eig2)
        log_eig2.sort()
        return log_eig2[::-1][:max_states]
    else:
        raise NotImplementedError("Only d=1 and d=2 are implemented for numerical eigenvalues.")

class NumericalEigenvalueLogW:
    def __init__(self, tau, N, d, L=2.0, num_grid=800, max_states=400):
        self.log_eigvals = get_nd_log_eigenvalues(tau, N, d, L, num_grid, max_states)
        self.N = N
        self.max_shell = len(self.log_eigvals) - 1
        
    def gk_fn(self, k):
        return 1
        
    def logwk_fn(self, k, logb):
        return self.log_eigvals[k] * self.N

def fermion_logZ_numeric_pimc_predict(tau, N, n, d=None, **kwargs):
    # Calculate eigenvalues dynamically for the given tau and N
    num_eig = NumericalEigenvalueLogW(tau, N, d, L=L_box, num_grid=800, max_states=max(n*4, 100))
    
    # Remove kwargs that conflict
    kwargs.pop('k_start', None)
    kwargs.pop('max_shell', None)
    kwargs.pop('gk_fn', None)
    kwargs.pop('logwk_fn', None)
    kwargs.pop('Ek_fn', None)
    
    return fermion_logZ_numeric(
        tau=tau, N=N, n=n, d=d, k_start=0, max_shell=num_eig.max_shell,
        gk_fn=num_eig.gk_fn, logwk_fn=num_eig.logwk_fn,
        **kwargs
    )
"""

# Insert eigen_code cell after Ek_box cell
new_cells = []
for cell in nb['cells']:
    new_cells.append(cell)
    if cell['cell_type'] == 'code' and 'Ek_box' in "".join(cell['source']) and 'val =' not in "".join(cell['source']):
        # Actually in my previous patch, Ek_box was defined in the same cell as val =.
        pass

# Let's just find the cell with Ek_box and val =
for cell in new_cells:
    if cell['cell_type'] == 'code':
        source = "".join(cell['source'])
        if "val =" in source and "Ek_box" in source:
            new_source = """n = 5

L_box = 2.0
def Ek_box(k):
    return (math.pi**2 / (2 * L_box**2)) * k

""" + eigen_code + """

val = (
    fermion_logZ_numeric_pimc_predict(tau=2.0, N=4, n=n, d=2)
    - fermion_logZ_numeric_pimc_predict(tau=2.5, N=4, n=n, d=2)
) / 0.5

print(val)
"""
            lines = [line + "\n" for line in new_source.split('\n')]
            if lines: lines[-1] = lines[-1][:-1]
            cell['source'] = lines
            cell['outputs'] = []

        elif "def plot_logZ_and_fd_multiN(" in source:
            # Add logZ_func to signature
            new_source = source.replace(
                "show_logZ=True, show_fd=True, **logZ_kwargs):",
                "show_logZ=True, show_fd=True, logZ_func=fermion_logZ_numeric, **logZ_kwargs):"
            ).replace(
                "fermion_logZ_numeric(tau=tau, N=N, n=n, d=d, **logZ_kwargs)",
                "logZ_func(tau=tau, N=N, n=n, d=d, **logZ_kwargs)"
            )
            cell['source'] = new_source.splitlines(keepends=True)

        elif "results = plot_logZ_and_fd_multiN(" in source:
            new_source = """n = 5
results = plot_logZ_and_fd_multiN(
    tau_start=1.0,
    tau_end=15.0,
    tau_step=1.0,
    
    n=n,
    d=2,
    N_list=[4],
    show_logZ=False,
    show_fd=True,
    logZ_func=fermion_logZ_numeric_pimc_predict
)
"""
            lines = [line + "\n" for line in new_source.split('\n')]
            if lines: lines[-1] = lines[-1][:-1]
            cell['source'] = lines
            cell['outputs'] = []

with open('JonasPolyConv-Copy3.ipynb', 'w') as f:
    json.dump(nb, f, indent=1)

