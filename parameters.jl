# --- Simulation Parameters ---
num_fermions = 400
dimensions = 2
bead_counts = [2, 4, 8, 16]
tau_start = 0.25
tau_stop = 15.25
tau_values = range(tau_start, tau_stop, length=60)

# --- High-Precision Settings ---
# Set precision for BigFloat (in bits). 50 decimal digits ≈ 167 bits.
bigfloat_precision = 170000

#Propagator choice: PA, TI, 4A is written as FA
propagator_choice = Symbol("PA")

mode = "ham" # thermo, ham

w = BigFloat(1/1) #For no interaction, set w=1
balanced = false

run_id = "comp_$mode _N$num_fermions _D$dimensions _$propagator_choice"

#Preview plot
toplim = 20000
bottomlim = 5000