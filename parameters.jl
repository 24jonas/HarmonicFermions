# --- Simulation Parameters ---
num_fermions = 3
dimensions = 2
bead_counts = [2, 4, 8, 12, 16]

tau_start = 0.25
tau_stop = 140.25
tau_values = range(tau_start, tau_stop, length=200)

# --- High-Precision Settings ---
# Set precision for BigFloat (in bits). 50 decimal digits ≈ 167 bits.
bigfloat_precision = 400

#Propagator choice: PA, TI, 4A is written as FA
propagator_choice = Symbol("FA")

mode = "ham" # thermo, ham

w = BigFloat(1/2) #For no interaction, set w=1
balanced = false

run_id = "normal"

#Preview plot
toplim =25
bottomlim = 0