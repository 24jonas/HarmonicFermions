using ArbNumerics

# --- Simulation Parameters ---
num_fermions = 100 #! cost
dimensions = 2
bead_counts = [2, 4, 8, 16] #! cost
tau_start = 0.5
tau_stop = 10.5 #! cost
tau_values = range(tau_start, tau_stop, length = 20) #! resolution (proportional cost)

# --- High-Precision Settings ---
bigfloat_precision = 20000 #! cost, determined by n, tau and to lesser extent N
setprecision(ArbFloat, bigfloat_precision)

#Propagator choice: PA, TI, 4A is written as FA
propagator_choice = Symbol("PA")

# Use ArbFloat explicitly to prevent type promotion errors later
w = ArbFloat(1.0) # Harmonic interactions parameter
balanced = false

# --- Output Directory ---
output_dir = "Output"
mkpath(output_dir)

run_id_T = "comp_thermo _N$num_fermions _D$dimensions _$propagator_choice"
run_id_H = "comp_ham _N$num_fermions _D$dimensions _$propagator_choice"
