# --- Simulation Parameters ---
num_fermions = 4
dimensions = 1
bead_counts = [2, 4, 6, 8, 200]

tau_start = 0.25
tau_stop = 10.25
tau_length = tau_stop - tau_start 

tau_values = range(tau_start, tau_stop, length=100)
toplim = 15
bottomlim = 0

# --- High-Precision Settings ---
# Set precision for BigFloat (in bits). 50 decimal digits ≈ 167 bits.
bigfloat_precision =  2000 #4096

#Propagator choice: PA, TI, 4A is written as FA
propagator_choice = Symbol("PA")

mode = "thermo" # thermo, ham

run_id = "normal"