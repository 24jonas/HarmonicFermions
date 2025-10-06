# --- Simulation Parameters ---
num_fermions = 6
dimensions = 2
bead_counts = [12, 27, 200]

tau_start = 5.25
tau_stop = 25.25
tau_values = range(tau_start, tau_stop, length=100)

# --- High-Precision Settings ---
# Set precision for BigFloat (in bits). 50 decimal digits ≈ 167 bits.
bigfloat_precision =  6000 #4096

#Propagator choice: PA, TI, 4A is written as FA
propagator_choice = Symbol("PA")

mode = "ham" # thermo, ham

w = BigFloat(1/2) #For no interaction, set w=1
balanced = true

run_id = "normal"

#Preview plot
toplim = 31
bottomlim = 5