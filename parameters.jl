# --- Simulation Parameters ---
num_fermions = 100
dimensions = 2
bead_counts = [2, 4, 8, 16]
tau_start = 0.25
tau_stop = 20.25
tau_values = range(tau_start, tau_stop, length=40)

# --- High-Precision Settings ---
# Set precision for BigFloat (in bits). 50 decimal digits ≈ 167 bits.
bigfloat_precision = 25000

#Propagator choice: PA, TI, 4A is written as FA
propagator_choice = Symbol("FA")

mode = "thermo" # thermo, ham

w = BigFloat(1/1) #For no interaction, set w=1
balanced = false

run_id = "comp_$mode _N$num_fermions _D$dimensions _$propagator_choice"

#Preview plot
toplim = 1200
bottomlim = 350