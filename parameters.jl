# --- Simulation Parameters ---
num_fermions = 3
dimensions = 3
bead_counts = [256]

tau_start = 0.25
tau_stop = 20.25
tau_values = range(tau_start, tau_stop, length=200)

# --- High-Precision Settings ---
# Set precision for BigFloat (in bits). 50 decimal digits ≈ 167 bits.
bigfloat_precision = 512

#Propagator choice: PA, TI, 4A is written as FA
propagator_choice = Symbol("FA")

mode = "exact" # thermo, ham

w = BigFloat(1/1) #For no interaction, set w=1
balanced = false

run_id = "$mode _N$num_fermions D$dimensions"

#Preview plot
toplim =15
bottomlim = 0