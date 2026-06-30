using ArbNumerics

# --- Simulation Parameters ---
num_fermions = 300
dimensions = 2
bead_counts = [2, 4, 8, 16]
tau_start = 0.5
tau_stop = 10.5
tau_values = range(tau_start, tau_stop, length = 20)

# --- High-Precision Settings ---
bigfloat_precision = 256 #Requires more precision than float64 but settles to a low value (relative to Direct recursion which keeps growing).
setprecision(ArbFloat, bigfloat_precision) #The extra precision is not needed if log exp sum trick is used, at the cost of slightly higher memory.

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
