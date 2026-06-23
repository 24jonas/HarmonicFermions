# HarmonicFermions

Exact computation of canonical partition functions, energies, and thermodynamic properties for **non-interacting fermions in harmonic traps**, using two complementary analytical methods.

## Overview

This repository implements two methods for computing the canonical partition function $Z_n(\beta)$ of $n$ non-interacting fermions:

1. **Baxter Method** — A combinatorial recursion based on Baxter's $q$-series identity, specialized for 2D harmonic oscillators. Computes $\log Q_n(b)$ where $b = e^{-\beta\omega}$ is the Boltzmann factor.

2. **Polynomial Convolution (PolyConv) Method** — Convolves per-shell generating functions $(1 + w_k t)^{g_k}$ in log-space to compute $\log Z_n$. Works for arbitrary dimensions and arbitrary single-particle potentials (not restricted to harmonic traps).

Both methods compute **exact** partition functions (up to floating-point precision) and are designed for use in **Path Integral Monte Carlo (PIMC)** simulations, where the ratio of fermionic to bosonic partition functions determines the sign of permutation weights.

## Repository Structure

### Shared Python Modules

| Module | Description |
|---|---|
| `baxter_method.py` | Core Baxter recursion: `log_Q(n, b)` and `log_Q_all(n, b)` with SQLite checkpointing |
| `polyconv_method.py` | General PolyConv: accepts custom degeneracy/weight functions for arbitrary potentials |
| `polyconv_opt.py` | Numba-accelerated PolyConv for harmonic spectra with SQLite checkpointing |
| `propagator.py` | Primitive Approximation (PA) propagator functions and energy estimator correction factors |
| `physics.py` | Exact ground-state energies (`energy_2d`, `energy_nd`), Thomas-Fermi approximation |
| `plotting.py` | Utilities: tau grids, finite differences, heat capacity computation |
| `benchmarking.py` | Runtime scaling analysis vs particle number, temperature, and bead count |

### Jupyter Notebooks

#### Baxter Method

| Notebook | Purpose |
|---|---|
| `BaxterMeth-opt.ipynb` | Energy estimators (thermodynamic & Hamiltonian) vs $\tau$ for 2D harmonic traps |
| `BaxterMeth-opt_freeTaus.ipynb` | Same as above with user-specified (non-uniform) $\tau$ grid |
| `BaxterMeth-heatmap.ipynb` | Heatmap: minimum beads $N$ for convergence of chemical potential |
| `BaxterMeth-heatmap-energy.ipynb` | Heatmap: minimum beads $N$ for convergence of energy |
| `BaxterMeth-heatmap-heatcap.ipynb` | Heatmap: minimum beads $N$ for convergence of heat capacity |

#### Polynomial Convolution Method

| Notebook | Purpose |
|---|---|
| `JonasPolyConv-opt.ipynb` | Numba-optimized PolyConv for harmonic traps with energy estimators |
| `JonasPolyConv-ArbitraryPotential.ipynb` | PolyConv for user-defined separable potentials (1D/2D) |
| `JonasPolyConv-SphericalPotential.ipynb` | PolyConv for user-defined spherically symmetric potentials (3D) |

#### Analysis

| Notebook | Purpose |
|---|---|
| `E_related_properties.ipynb` | Energy, heat capacity, and chemical potential vs $n$; Thomas-Fermi comparison |

### Other

| Path | Description |
|---|---|
| `Saved_runs_and_plots/` | Pre-computed data (CSV, EPS) and generated plots |
| `Julia_PolyConv/` | Julia implementation of the PolyConv method (QMC simulation) |

## Getting Started

### Requirements

```
numpy
numba
scipy
matplotlib
pandas
joblib
```

### Installation

```bash
git clone https://github.com/<your-username>/HarmonicFermions.git
cd HarmonicFermions
pip install numpy numba scipy matplotlib pandas joblib
```

### Quick Start

**Compute the partition function for 100 fermions in a 2D harmonic trap:**

```python
from polyconv_opt import fermion_logZ_numeric

# tau = inverse temperature * bead count, N = number of Trotter beads
logZ = fermion_logZ_numeric(tau=5.0, N=4, n=100, d=2, w=1.0)
print(f"log Z_100 = {logZ:.6f}")
```

**Compute the exact ground-state energy:**

```python
from physics import energy_2d, energy_nd

E_2d = energy_2d(n=100, w=1.0)      # 2D harmonic trap
E_3d = energy_nd(n=100, d=3, w=1.0)  # 3D harmonic trap
print(f"E(n=100, 2D) = {E_2d}")
print(f"E(n=100, 3D) = {E_3d}")
```

**Use the Baxter method:**

```python
from baxter_method import log_Q
from propagator import get_b_val

b = get_b_val(tau_val=5.0, N_val=4, w_val=1.0)
logQ = log_Q(n=100, b=b)
print(f"log Q_100 = {logQ:.6f}")
```

**Compute partition functions with a custom potential:**

```python
import numpy as np
from polyconv_method import fermion_logZ_numeric

# Define a 1D double-well potential (requires numerical eigenvalue setup,
# see JonasPolyConv-ArbitraryPotential.ipynb for full examples)
logZ = fermion_logZ_numeric(tau=2.0, N=4, n=10, d=2)
```

For full examples with custom potentials, see the `JonasPolyConv-ArbitraryPotential.ipynb` and `JonasPolyConv-SphericalPotential.ipynb` notebooks.

## Methods

### Baxter Method

Uses the identity relating the canonical partition function of 2D fermions to Baxter's $q$-series. The recursion computes in $O(n^2)$ time with $O(n)$ memory, and supports checkpoint/resume via SQLite for large $n$.

### Polynomial Convolution

Builds the $n$-particle partition function by multiplying per-shell generating functions:

$$Z_n = [t^n] \prod_{k=0}^{\infty} (1 + w_k\, t)^{g_k}$$

where $g_k$ is the degeneracy and $w_k = e^{-\beta E_k}$ the Boltzmann weight of shell $k$. All arithmetic is done in log-space for numerical stability. The computation adapts to the spectrum and converges once the contributions of higher shells become negligible.

### Energy Estimators

Two energy estimators are computed from $\log Z_n(\tau)$:

- **Thermodynamic estimator**: $E_T = -\partial \log Z / \partial \tau$ (via finite difference)
- **Hamiltonian estimator**: $E_H = E_T \times f_H / f_T$ (corrected by propagator-dependent factors)

Both converge to the exact energy as $N \to \infty$ (Trotter limit).

## License

This project is provided as-is for academic use.
