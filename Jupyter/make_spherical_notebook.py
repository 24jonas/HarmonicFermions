import json
import shutil

shutil.copyfile('JonasPolyConv-ArbitraryPotential.ipynb', 'JonasPolyConv-SphericalPotential.ipynb')

with open('JonasPolyConv-SphericalPotential.ipynb', 'r') as f:
    nb = json.load(f)

# The user input cell is cell 1
# Let's replace its content
user_input_content = [
    "import numpy as np\n",
    "\n",
    "# ==========================================\n",
    "# USER INPUTS: Define your spherical potential\n",
    "# ==========================================\n",
    "\n",
    "# Radial Integration Domain\n",
    "r_min = 1e-6\n",
    "r_max = 5.0\n",
    "num_grid = 800\n",
    "\n",
    "# Maximum angular momentum to compute\n",
    "l_max = 10\n",
    "\n",
    "def user_potential_radial(r):\n",
    "    \"\"\"\n",
    "    Define the spherically symmetric potential V(r).\n",
    "    \"\"\"\n",
    "    # --- Default: Regularized Hydrogen Atom ---\n",
    "    # V(r) = -1 / sqrt(r^2 + a^2)\n",
    "    a = 0.5\n",
    "    return -1.0 / np.sqrt(r**2 + a**2)\n",
    "\n",
    "    # --- Example 2: 3D Isotropic Harmonic Oscillator ---\n",
    "    # w = 1.0\n",
    "    # return 0.5 * (w**2) * (r**2)\n",
    "\n",
    "    # --- Example 3: 3D Spherical Well ---\n",
    "    # V = np.zeros_like(r)\n",
    "    # V[r > 2.0] = np.inf\n",
    "    # return V\n",
    "\n",
    "# Parameters for testing\n",
    "n = 5\n",
    "N = 4\n"
]

nb['cells'][1]['source'] = user_input_content

# The solver code is in cell 2
# Let's replace it
solver_content = """import numpy as np
import scipy.linalg as la
import scipy.special

def get_radial_log_eigenvalues(tau, N, r_min, r_max, num_grid, l_max):
    epsilon = tau / N
    r = np.linspace(r_min, r_max, num_grid)
    dr = (r_max - r_min) / (num_grid - 1)
    
    R, Rp = np.meshgrid(r, r)
    z = R * Rp / epsilon
    
    V_R = user_potential_radial(R)
    V_Rp = user_potential_radial(Rp)
    
    all_log_eigvals = []
    
    # Calculate for each angular momentum l
    for l in range(l_max + 1):
        # Exact 3D radial primitive approximation density matrix
        rho_l = (R * Rp / epsilon) * np.exp(-((R - Rp)**2) / (2 * epsilon)) * scipy.special.ive(l + 0.5, z)
        T_l = np.exp(-epsilon * V_R / 2.0) * rho_l * np.exp(-epsilon * V_Rp / 2.0) * dr
        
        eigvals = la.eigvalsh(T_l)
        eigvals = eigvals[eigvals > 1e-15] # filter numerical noise
        log_eigvals_l = np.log(eigvals)
        
        # Each eigenvalue for angular momentum l has degeneracy 2l + 1
        for log_val in log_eigvals_l:
            all_log_eigvals.append((log_val, 2 * l + 1))
            
    # Sort all states descending by eigenvalue
    all_log_eigvals.sort(key=lambda x: x[0], reverse=True)
    return all_log_eigvals

class NumericalEigenvalueLogW:
    def __init__(self, tau, N, r_min, r_max, num_grid, l_max):
        # list of (log_eigval, degeneracy) tuples
        self.sorted_states = get_radial_log_eigenvalues(tau, N, r_min, r_max, num_grid, l_max)
        self.N = N
        self.max_shell = len(self.sorted_states) - 1
        
    def gk_fn(self, k):
        return self.sorted_states[k][1]
        
    def logwk_fn(self, k, logb):
        return self.sorted_states[k][0] * self.N

def fermion_logZ_numeric_pimc_predict(tau, N, n, d=None, **kwargs):
    num_eig = NumericalEigenvalueLogW(tau, N, r_min, r_max, num_grid, l_max)
    
    kwargs.pop('k_start', None)
    kwargs.pop('max_shell', None)
    kwargs.pop('gk_fn', None)
    kwargs.pop('logwk_fn', None)
    kwargs.pop('Ek_fn', None)
    
    return fermion_logZ_numeric(
        tau=tau, N=N, n=n, d=None, k_start=0, max_shell=num_eig.max_shell,
        gk_fn=num_eig.gk_fn, logwk_fn=num_eig.logwk_fn,
        **kwargs
    )
"""

for i, cell in enumerate(nb['cells']):
    if cell['cell_type'] == 'code':
        source = "".join(cell['source'])
        if "def get_1d_log_eigenvalues" in source:
            lines = [line + "\n" for line in solver_content.split('\n')]
            if lines: lines[-1] = lines[-1][:-1]
            cell['source'] = lines

        elif "val =" in source and "fermion_logZ_numeric_pimc_predict" in source:
            new_source = """val = (
    fermion_logZ_numeric_pimc_predict(tau=2.0, N=4, n=n)
    - fermion_logZ_numeric_pimc_predict(tau=2.5, N=4, n=n)
) / 0.5

print(val)
"""
            lines = [line + "\n" for line in new_source.split('\n')]
            if lines: lines[-1] = lines[-1][:-1]
            cell['source'] = lines

        elif "plot_logZ_and_fd_multiN(" in source and "d=2" in source:
            new_source = source.replace("d=2", "d=None")
            cell['source'] = new_source.splitlines(keepends=True)


with open('JonasPolyConv-SphericalPotential.ipynb', 'w') as f:
    json.dump(nb, f, indent=1)

