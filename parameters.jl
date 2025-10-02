# --- Simulation Parameters ---
num_fermions = 6
dimensions = 2
bead_counts = [12, 18, 27, 40, 200]

tau_start = 5.25
tau_stop = 30.25
tau_length = tau_stop - tau_start 

tau_values = range(tau_start, tau_stop, length=100)
toplim = 5.625
bottomlim = 5.475

# --- High-Precision Settings ---
# Set precision for BigFloat (in bits). 50 decimal digits ≈ 167 bits.
bigfloat_precision =  2000 #4096

#Propagator choice: PA, TI, 4A is written as FA
propagator_choice = Symbol("PA")

mode = "ham" # thermo, ham

run_id = "normal"