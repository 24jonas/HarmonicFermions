# ----------------------------------------------------
# Parameters
# ----------------------------------------------------
n = 4        # Number of particles
N = 8       # Number of beads

# ----------------------------------------------------
# EPS Export Settings
# ----------------------------------------------------
# Set terminal to EPS, enhanced text, color, and size
set terminal postscript eps enhanced color font "Helvetica,14" linewidth 4 size 6,4.5
set output "energy_heat_capacity_overlay.eps"

# ----------------------------------------------------
# Helper Functions
# ----------------------------------------------------
# Hyperbolic functions
coth(x) = cosh(x) / sinh(x)
csch(x) = 1.0 / sinh(x)

# Updated exact helper functions for H and T variants
Ex(tau) = tau / N
g(tau)  = sqrt(1.0 + (Ex(tau)**2) / 4.0)
Z(tau)  = 1.0 + (Ex(tau)**2) / 2.0
u(tau)  = acosh(Z(tau))
b(tau)  = exp(-N * u(tau))

# ----------------------------------------------------
# Exact Energy and Heat Capacity Equations
# ----------------------------------------------------
# Exact Energy E_n(x)
E_n(tau) = n*(n-1)/4.0 + 0.5 * sum [i=1:n] (i * coth(i * N * asinh(tau/(2.0*N))))

# Exact Heat Capacity C_n(x)
C_n_exact(tau) = sum [i=1:n] ( (i**2 * tau**2 / 4.0) * (csch(i*tau/2.0)**2) )

# ----------------------------------------------------
# Approximate Energies
# ----------------------------------------------------
# E_TnA(x) (using 1.0 for the 'l' in the prompt's numerator)
E_TnA(tau) = (1.0 / g(tau)) * ( n**2 / 2.0 + sum [i=1:n] ( (i * b(tau)**i) / (1.0 - b(tau)**i) ) )

# E_HnA(x)
E_HnA(tau) = 0.5 * (1.0/g(tau) + g(tau)) * ( n**2 / 2.0 + sum [i=1:n] ( (i * b(tau)**i) / (1.0 - b(tau)**i) ) )

# ----------------------------------------------------
# Numerical Derivatives for Approximated Heat Capacities
# ----------------------------------------------------
delta = 1e-6 

dE_TnA_dtau(tau) = (E_TnA(tau + delta) - E_TnA(tau - delta)) / (2.0 * delta)
dE_HnA_dtau(tau) = (E_HnA(tau + delta) - E_HnA(tau - delta)) / (2.0 * delta)

# Heat capacities C(x) = -x^2 * (dE/dx) / n
C_TnA(tau) = -tau**2 * dE_TnA_dtau(tau) / n
C_HnA(tau) = -tau**2 * dE_HnA_dtau(tau) / n

# ----------------------------------------------------
# Plot Settings
# ----------------------------------------------------
# set title "Thermodynamics for n=".n.", N=".N font "Helvetica,16"

# Common styling and grid
set xlabel "T" font "Helvetica,20"
set ylabel "Heat Capacity / n   &   2 * Energy / n^2" font "Helvetica,20"
set grid linetype 1 linecolor rgb "#cccccc" linewidth 1
set samples 1000

# Hide default legend since we are using explicit labels
unset key 

# T domain
T_min = 0.01 
T_max = 2.0

set yrange [0:1.5]
# Added font sizing to tics to make the grid numbers larger
#set xtics 0.25 font "Helvetica,18"
#set ytics 0.25 font "Helvetica,18"

# ----------------------------------------------------
# Custom Labels (Adjust coordinates as needed)
# ----------------------------------------------------
# Energy Labels
set label "E"   at 1.7, 1.25 font "Helvetica,20" textcolor rgb "black"
set label "E_T" at 0.6, 0.95 font "Helvetica,20" textcolor rgb "red"
set label "E_H" at 0.6, 1.15 font "Helvetica,20" textcolor rgb "blue"

# Heat Capacity Labels
set label "C"   at 1.7, 0.7 font "Helvetica,20" textcolor rgb "black"
set label "C_T" at 0.6, 0.5 font "Helvetica,20" textcolor rgb "red"
set label "C_H" at 0.6, 0.25 font "Helvetica,20" textcolor rgb "blue"

# ----------------------------------------------------
# Overlay Plot Command
# ----------------------------------------------------
plot [T=T_min:T_max] \
     E_n(1.0/T) / (0.5*n**2)     title "Exact E"   with lines lw 3 lc rgb "black", \
     E_TnA(1.0/T) / (0.5*n**2)   title "E_T"       with lines lw 3 dt 2 lc rgb "red", \
     E_HnA(1.0/T) / (0.5*n**2)   title "E_H"       with lines lw 3 dt 3 lc rgb "blue", \
     C_n_exact(1.0/T) / n     title "Exact C"   with lines lw 3 lc rgb "black", \
     C_TnA(1.0/T)             title "C_T"       with lines lw 3 dt 2 lc rgb "red", \
     C_HnA(1.0/T)             title "C_H"       with lines lw 3 dt 3 lc rgb "blue"
