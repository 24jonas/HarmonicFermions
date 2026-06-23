# ----------------------------------------------------
# Parameters
# ----------------------------------------------------
n = 50        # Number of particles

# ----------------------------------------------------
# EPS Export Settings (Common)
# ----------------------------------------------------
set terminal postscript eps enhanced color font "Helvetica,14" linewidth 4 size 6,4.5

# ----------------------------------------------------
# Helper Functions
# ----------------------------------------------------
coth(x) = cosh(x) / sinh(x)
csch(x) = 1.0 / sinh(x)

Ex(tau, N) = tau / N

# --- PA Propagator Helpers ---
g(tau, N)  = sqrt(1.0 + (Ex(tau, N)**2) / 4.0)
Z(tau, N)  = 1.0 + (Ex(tau, N)**2) / 2.0
u(tau, N)  = acosh(Z(tau, N))
b(tau, N)  = exp(-N * u(tau, N))

# --- 4A (FA) Propagator Helpers ---
Z_4A(tau, N) = 1.0 + (Ex(tau, N)**2)/2.0 + (Ex(tau, N)**4)/24.0 + (Ex(tau, N)**6)/864.0
g_4A(tau, N) = sqrt(1.0 + (Ex(tau, N)**4) / (432.0 + 36.0*(Ex(tau, N)**2)))
u_4A(tau, N) = acosh(Z_4A(tau, N))
b_4A(tau, N) = exp(-N * u_4A(tau, N))

# ----------------------------------------------------
# Exact Energy and Heat Capacity Equations
# ----------------------------------------------------
E_n(tau, N) = n*(n-1)/4.0 + 0.5 * sum [i=1:n] (i * coth(i * N * asinh(tau/(2.0*N))))
C_n_exact(tau) = sum [i=1:n] ( (i**2 * tau**2 / 4.0) * (csch(i*tau/2.0)**2) )

# ----------------------------------------------------
# Approximate Energies
# ----------------------------------------------------
# PA Approximate Energies
E_TnA(tau, N) = (1.0 / g(tau, N)) * ( n**2 / 2.0 + sum [i=1:n] ( (i * b(tau, N)**i) / (1.0 - b(tau, N)**i) ) )
E_HnA(tau, N) = 0.5 * (1.0/g(tau, N) + g(tau, N)) * ( n**2 / 2.0 + sum [i=1:n] ( (i * b(tau, N)**i) / (1.0 - b(tau, N)**i) ) )

# 4A Approximate Energies
E_TnA_4A(tau, N) = (1.0 / g_4A(tau, N)) * ( n**2 / 2.0 + sum [i=1:n] ( (i * b_4A(tau, N)**i) / (1.0 - b_4A(tau, N)**i) ) )
E_HnA_4A(tau, N) = 0.5 * (1.0/g_4A(tau, N) + g_4A(tau, N)) * ( n**2 / 2.0 + sum [i=1:n] ( (i * b_4A(tau, N)**i) / (1.0 - b_4A(tau, N)**i) ) )

# ----------------------------------------------------
# Numerical Derivatives for Approximated Heat Capacities
# ----------------------------------------------------
delta = 1e-6

# PA Heat Capacities
dE_TnA_dtau(tau, N) = (E_TnA(tau + delta, N) - E_TnA(tau - delta, N)) / (2.0 * delta)
dE_HnA_dtau(tau, N) = (E_HnA(tau + delta, N) - E_HnA(tau - delta, N)) / (2.0 * delta)
C_TnA(tau, N) = -tau**2 * dE_TnA_dtau(tau, N) / n
C_HnA(tau, N) = -tau**2 * dE_HnA_dtau(tau, N) / n

# 4A Heat Capacities
dE_TnA_4A_dtau(tau, N) = (E_TnA_4A(tau + delta, N) - E_TnA_4A(tau - delta, N)) / (2.0 * delta)
dE_HnA_4A_dtau(tau, N) = (E_HnA_4A(tau + delta, N) - E_HnA_4A(tau - delta, N)) / (2.0 * delta)
C_TnA_4A(tau, N) = -tau**2 * dE_TnA_4A_dtau(tau, N) / n
C_HnA_4A(tau, N) = -tau**2 * dE_HnA_4A_dtau(tau, N) / n

# ----------------------------------------------------
# Common Plot Settings
# ----------------------------------------------------
set xlabel "{{/Symbol t}" font "Helvetica,20"
set grid linetype 1 linecolor rgb "#cccccc" linewidth 1
set samples 1000
set title "n=50, N=[4, 8, 16]" font "Helvetica,24"
unset key

T_min = 0.01
T_max = 15.0

# ----------------------------------------------------
# PLOT 1: ENERGIES ONLY
# ----------------------------------------------------
set output "energies_multi_N_lown.eps"
set ylabel "Energy" font "Helvetica,20"
set yrange [1250*0.999:1250*1.001]

set label "E_{continuum}"    at 13.1, 1249.9 font "Helvetica,18" textcolor rgb "black"
set label "T^{PA}"      at 2, 1249.3 font "Helvetica,18" textcolor rgb "blue"
set label "H^{PA}"      at 5, 1250.6 font "Helvetica,18" textcolor rgb "red"
set label "T^{4A}" at 7, 1249.6 font "Helvetica,18" textcolor rgb "#9D00FF"
set label "H^{4A}" at 13, 1250.2 font "Helvetica,18" textcolor rgb "#FF6E0D"

plot [T=T_min:T_max] \
     E_TnA(T, 4)            with lines lw 2  lc rgb "blue", \
     E_TnA(T, 8)            with lines lw 2  lc rgb "blue", \
     E_TnA(T, 16)          with lines lw 2  lc rgb "blue", \
     E_HnA(T, 4)            with lines lw 2  lc rgb "red", \
     E_HnA(T, 8)            with lines lw 2  lc rgb "red", \
     E_HnA(T, 16)          with lines lw 2  lc rgb "red", \
     E_TnA_4A(T, 4)         with lines lw 2  lc rgb "#9D00FF", \
     E_TnA_4A(T, 8)         with lines lw 2  lc rgb "#9D00FF", \
     E_TnA_4A(T, 16)       with lines lw 2  lc rgb "#9D00FF", \
     E_HnA_4A(T, 4)         with lines lw 2  lc rgb "#FF6E0D", \
     E_HnA_4A(T, 8)         with lines lw 2  lc rgb "#FF6E0D", \
     E_HnA_4A(T, 16)       with lines lw 2  lc rgb "#FF6E0D", \
     E_n(T, 32)           with lines lw 2 dt 2 lc rgb "black"

# ----------------------------------------------------
# PLOT 2: HEAT CAPACITIES ONLY
# ----------------------------------------------------
set output "heat_capacities_multi_N_lown.eps"
set ylabel "Heat Capacity / n" font "Helvetica,20"
set yrange [-0.5:1.5]
set xlabel "T"


unset label
set label "C_{continuum}"    at 1.7, 0.03 font "Helvetica,18" textcolor rgb "black"
set label "T^{PA}"      at 0.85, 0.75 font "Helvetica,18" textcolor rgb "blue"
set label "H^{PA}"      at 0.45, -0.2 font "Helvetica,18" textcolor rgb "red"
set label "T^{4A}" at 0.3, 0.30 font "Helvetica,18" textcolor rgb "#9D00FF"
set label "H^{4A}" at 0.01, 0.05 font "Helvetica,18" textcolor rgb "#FF6E0D"

plot [T=T_min:2] \
     C_TnA(1.0/T, 4)                with lines lw 2  lc rgb "blue", \
     C_TnA(1.0/T, 8)                with lines lw 2  lc rgb "blue", \
     C_TnA(1.0/T, 16)               with lines lw 2  lc rgb "blue", \
     C_HnA(1.0/T, 4)                 with lines lw 2  lc rgb "red", \
     C_HnA(1.0/T, 8)                 with lines lw 2  lc rgb "red", \
     C_HnA(1.0/T, 16)               with lines lw 2  lc rgb "red", \
     C_TnA_4A(1.0/T, 4)            with lines lw 2  lc rgb "#9D00FF", \
     C_TnA_4A(1.0/T, 8)             with lines lw 2  lc rgb "#9D00FF", \
     C_TnA_4A(1.0/T, 16)            with lines lw 2  lc rgb "#9D00FF", \
     C_HnA_4A(1.0/T, 4)              with lines lw 2  lc rgb "#FF6E0D", \
     C_HnA_4A(1.0/T, 8)              with lines lw 2  lc rgb "#FF6E0D", \
     C_HnA_4A(1.0/T, 16)            with lines lw 2  lc rgb "#FF6E0D", \
     C_n_exact(1.0/T) / n         with lines lw 2 dt 2 lc rgb "black"
