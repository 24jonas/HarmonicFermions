# HarmonicFermions

Exact computation of partition functions, energies, and thermodynamic properties for **non-interacting fermions in harmonic traps**, using three complementary methods.

## Overview

This repository implements three methods for computing the canonical partition function $Z_n(\beta)$ of $n$ non-interacting fermions:

1. **Direct Recursion** (Julia) — The Newton–Girard recursion for $Z_n$ from single-particle partition functions $z_k$, using arbitrary-precision arithmetic (`ArbNumerics.jl`). Supports multiple higher-order propagators (PA, TI, 4th-order) and computes both thermodynamic and Hamiltonian energy estimators.

2. **Baxter Method** (Python) — A combinatorial recursion based on Baxter's $q$-series identity, specialized for 2D harmonic oscillators. Computes $\log Q_n(b)$ where $b = e^{-\beta\omega}$ is the Boltzmann factor.

3. **Polynomial Convolution (PolyConv)** — Convolves per-shell generating functions $(1 + w_k t)^{g_k}$ in log-space to compute $\log Z_n$. Works for arbitrary dimensions and arbitrary single-particle potentials (not restricted to harmonic traps). Available in both Python (with Numba acceleration) and Julia (with `ArbNumerics` precision). *In development

All three methods compute **exact** partition functions (up to floating-point or arbitrary precision) and are designed for use in **Path Integral Monte Carlo (PIMC)** simulations, where the ratio of fermionic to bosonic partition functions determines the sign of permutation weights.

## Repository Structure

```
HarmonicFermions/
├── Direct_recursion/            # Julia: Newton–Girard recursion (arbitrary precision)
│   ├── parameters_ARB.jl        # Simulation parameters and output directory
│   ├── qmc_simulation_ARB.jl    # Main simulation driver
│   ├── EnergySol_ARB.jl         # Energy computation via Newton–Girard recursion
│   ├── PropagatorsModule_ARB.jl # PA, TI, and 4th-order propagators
│   ├── heat_capacity.jl         # Heat capacity from finite differences
│   └── Plots and previous runs/ # Archived run data and plots
│
├── Baxter_and_PolyConv/         # Python & Julia: Baxter and PolyConv methods
│   ├── baxter_method.py         # Baxter q-series recursion with SQLite checkpointing
│   ├── polyconv_method.py       # General PolyConv for arbitrary potentials
│   ├── polyconv_opt.py          # Numba-accelerated PolyConv for harmonic spectra
│   ├── propagator.py            # PA propagator functions and correction factors
│   ├── physics.py               # Exact ground-state energies, Thomas-Fermi approx.
│   ├── plotting.py              # Tau grids, finite differences, heat capacity
│   ├── benchmarking.py          # Runtime scaling analysis (vs n, tau, N)
│   │
│   ├── BaxterMath-opt.ipynb              # Baxter: energy vs tau for 2D traps
│   ├── BaxterMeth-heatmap-combined.ipynb # Heatmaps: bead convergence for Energy, ChemPot, HeatCap
│   ├── PolyConv-opt.ipynb                # PolyConv: energy estimators (harmonic)
│   ├── PolyConv-ArbitraryPotential.ipynb # PolyConv: custom separable potentials
│   ├── PolyConv-SphericalPotential.ipynb # PolyConv: custom spherical potentials
│   ├── E_related_properties.ipynb        # Energy, heat capacity, chem. potential vs n
│   │
│   ├── Julia_PolyConv/          # Julia: PolyConv with arbitrary precision
│   │   ├── parameters_ARB.jl
│   │   ├── qmc_simulation_ARB.jl
│   │   ├── EnergySol_ARB.jl     # Grand-Z polynomial convolution approach
│   │   ├── PropagatorsModule_ARB.jl
│   │   └── heat_capacity.jl
│   │
│   └── Saved_runs_and_plots/    # Archived run data and plots
│
└── .gitignore
```

## Methods

### 1. Direct Recursion (Newton–Girard)

Computes $Z_n$ from single-particle partition functions using the Newton–Girard identity:

$$Z_n = \frac{1}{n} \sum_{k=1}^{n} (-1)^{k-1} z_k \, Z_{n-k}$$

where $z_k = \text{tr}(e^{-k\beta H_1})$ is the $k$-particle trace. The Julia implementation uses `ArbNumerics.jl` for arbitrary-precision arithmetic (typically 10,000–20,000 bits), which is necessary because the alternating-sign recursion causes severe cancellation at large $n$ and $\tau$.

**Supports three propagator approximations:** Primitive (PA), Takahashi–Imada (TI), and 4th-order (4A).

### 2. Baxter Method

Uses the identity relating the canonical partition function of 2D fermions to Baxter's $q$-series. The recursion computes in $O(n^2)$ time with $O(n)$ memory, and supports checkpoint/resume via SQLite for large $n$. Specialized to $d = 2$.

### 3. Polynomial Convolution (PolyConv)

Builds the $n$-particle partition function by multiplying per-shell generating functions:

$$Z_n = [t^n] \prod_{k=0}^{\infty} (1 + w_k\, t)^{g_k}$$

where $g_k$ is the degeneracy and $w_k = e^{-\beta E_k}$ the Boltzmann weight of shell $k$. All arithmetic is done in log-space for numerical stability. The computation adapts to the spectrum and converges once the contributions of higher shells become negligible.

The PolyConv method works in **arbitrary dimensions** and accepts **custom single-particle spectra**, making it applicable to non-harmonic potentials. The precision requirements are much lower than the direct recursion (typically 256 bits vs 20,000+).

### Energy Estimators

Two energy estimators are computed from $\log Z_n(\tau)$:

- **Thermodynamic estimator**: $E_T = -\partial \log Z / \partial \tau$ (via finite difference)
- **Hamiltonian estimator**: $E_H = E_T \times f_H / f_T$ (corrected by propagator-dependent factors)

Both converge to the exact energy as $N \to \infty$ (Trotter limit).

## Getting Started

### Python Requirements

```
numpy
numba
scipy
matplotlib
pandas
```

### Julia Requirements

```julia
using Pkg
Pkg.add(["ArbNumerics", "Polynomials", "DataFrames", "CSV", "JSON"])
```

### Quick Start

**Compute the partition function for 100 fermions in a 2D harmonic trap (Python):**

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

**Run the Direct Recursion (Julia):**

```bash
cd Direct_recursion
# Edit parameters_ARB.jl to set num_fermions, dimensions, bead_counts, etc.
julia qmc_simulation_ARB.jl
```

**Use the PolyConv method with a custom potential:**

```python
from polyconv_method import fermion_logZ_numeric

# See PolyConv-ArbitraryPotential.ipynb for full examples with
# custom degeneracy (gk_fn) and weight (logwk_fn) functions
logZ = fermion_logZ_numeric(tau=2.0, N=4, n=10, d=2)
```

### Output

All methods write results to an `Output/` directory (configurable via `output_dir` in parameters files and notebook cells). The Julia implementations use `CSV.jl`; the Python notebooks save `.csv` and `.npy` files.

## License

This project is provided as-is for academic use.
