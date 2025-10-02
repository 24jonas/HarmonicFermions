# --- Simulation Parameters ---
num_fermions = 6
dimensions = 2
bead_counts = [12, 27, 200]

tau_start = 0.25
tau_stop = 20.25
tau_length = tau_stop - tau_start 

tau_values = range(tau_start, tau_stop, length=100)
toplim = 5.62
bottomlim = 5.48

# --- High-Precision Settings ---
# Set precision for BigFloat (in bits). 50 decimal digits ≈ 167 bits.
bigfloat_precision =  2000 #4096

#Propagator choice: PA, TI, 4A is written as FA
propagator_choice = Symbol("PA")

mode = "ham" # thermo, ham

run_id = "normal"