# --- Simulation Parameters ---
num_fermions = 100
dimensions = 2
bead_counts = [200]
tau_start = 0.25
tau_stop = 10.25
tau_values = range(tau_start, tau_stop, length=200)

# --- High-Precision Settings ---
# Set precision for BigFloat (in bits). 50 decimal digits ≈ 167 bits.
bigfloat_precision = 20000

#Propagator choice: PA, TI, 4A is written as FA
propagator_choice = Symbol("PA")

mode = "exact" # thermo, ham

w = BigFloat(1/2) #For no interaction, set w=1
balanced = false

run_id = "comp_$mode _N$num_fermions _D$dimensions _$propagator_choice"

#Preview plot
toplim = 1000
bottomlim = 1