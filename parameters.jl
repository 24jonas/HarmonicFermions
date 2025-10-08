# --- Simulation Parameters ---
num_fermions = 3
dimensions = 2
bead_counts = [2, 4, 6, 8, 64]

tau_start = 0.25
tau_stop = 20.25
tau_values = range(tau_start, tau_stop, length=100)

# --- High-Precision Settings ---
# Set precision for BigFloat (in bits). 50 decimal digits ≈ 167 bits.
bigfloat_precision = 200

#Propagator choice: PA, TI, 4A is written as FA
propagator_choice = Symbol("PA")

mode = "ham" # thermo, ham

w = BigFloat(3/2) #For no interaction, set w=1
balanced = false

run_id = "normal"

#Preview plot
toplim = 12
bottomlim = 2