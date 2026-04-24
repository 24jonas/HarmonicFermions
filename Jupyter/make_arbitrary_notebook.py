import json
import shutil

# Copy notebook
shutil.copyfile('JonasPolyConv-Copy3.ipynb', 'JonasPolyConv-ArbitraryPotential.ipynb')

with open('JonasPolyConv-ArbitraryPotential.ipynb', 'r') as f:
    nb = json.load(f)

user_input_cell = {
    "cell_type": "code",
    "execution_count": None,
    "metadata": {},
    "outputs": [],
    "source": [
        "import numpy as np\n",
        "\n",
        "# ==========================================\n",
        "# USER INPUTS: Define your potential here\n",
        "# ==========================================\n",
        "\n",
        "# Integration Domain\n",
        "x_min = -1.0\n",
        "x_max = 1.0\n",
        "num_grid = 800\n",
        "\n",
        "def user_potential(x):\n",
        "    \"\"\"\n",
        "    Define the separable 1D potential V(x).\n",
        "    Must be able to accept a NumPy array x and return an array of the same shape.\n",
        "    \"\"\"\n",
        "    # --- Example 1: Hard Wall (Infinite Square Well) ---\n",
        "    # (V=0 inside domain, boundary is handled by x_min, x_max)\n",
        "    return np.zeros_like(x)\n",
        "\n",
        "    # --- Example 2: Harmonic Oscillator ---\n",
        "    # w = 1.0\n",
        "    # return 0.5 * (w**2) * (x**2)\n",
        "\n",
        "    # --- Example 3: Double Well ---\n",
        "    # return (x**2 - 1.0)**2\n",
        "\n",
        "# Parameters for testing\n",
        "n = 5\n",
        "d = 2\n",
        "N = 4\n"
    ]
}

# Find the cell containing get_1d_log_eigenvalues
for i, cell in enumerate(nb['cells']):
    if cell['cell_type'] == 'code':
        source = "".join(cell['source'])
        if "def get_1d_log_eigenvalues" in source:
            new_source = """import numpy as np
import scipy.linalg as la

def get_1d_log_eigenvalues(tau, N, x_min, x_max, num_grid):
    epsilon = tau / N
    x = np.linspace(x_min, x_max, num_grid)
    dx = (x_max - x_min) / (num_grid - 1)
    
    X, Y = np.meshgrid(x, x)
    
    # Calculate Potential Energy components
    Vx = user_potential(X)
    Vy = user_potential(Y)
    
    # Calculate Kinetic Matrix
    K_kinetic = (1.0 / np.sqrt(2 * np.pi * epsilon)) * np.exp(-((X - Y)**2) / (2 * epsilon))
    
    # Symmetric Trotter Breakup for Primitive Approximation
    T = np.exp(-epsilon * Vx / 2.0) * K_kinetic * np.exp(-epsilon * Vy / 2.0) * dx
    
    # Extract eigenvalues of the symmetric transfer matrix
    eigvals = la.eigvalsh(T)
    eigvals = eigvals[eigvals > 1e-15] # filter numerical noise
    return np.log(np.sort(eigvals)[::-1])

def get_nd_log_eigenvalues(tau, N, d, x_min, x_max, num_grid, max_states=1000):
    log_eig1 = get_1d_log_eigenvalues(tau, N, x_min, x_max, num_grid)
    
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
    def __init__(self, tau, N, d, x_min, x_max, num_grid, max_states=400):
        self.log_eigvals = get_nd_log_eigenvalues(tau, N, d, x_min, x_max, num_grid, max_states)
        self.N = N
        self.max_shell = len(self.log_eigvals) - 1
        
    def gk_fn(self, k):
        return 1
        
    def logwk_fn(self, k, logb):
        return self.log_eigvals[k] * self.N

def fermion_logZ_numeric_pimc_predict(tau, N, n, d=None, **kwargs):
    # Calculate eigenvalues dynamically for the given tau and N using user_potential
    num_eig = NumericalEigenvalueLogW(tau, N, d, x_min, x_max, num_grid, max_states=max(n*4, 100))
    
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
            # Also remove Ek_box definitions from this cell since we don't need it
            lines = [line + "\n" for line in new_source.split('\n')]
            if lines: lines[-1] = lines[-1][:-1]
            cell['source'] = lines
            break

# Prepend the user_input_cell at the top (after fermion_logZ_numeric definition, which is cell 0)
nb['cells'].insert(1, user_input_cell)

with open('JonasPolyConv-ArbitraryPotential.ipynb', 'w') as f:
    json.dump(nb, f, indent=1)

