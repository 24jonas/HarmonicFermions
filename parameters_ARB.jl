using ArbNumerics

# --- Simulation Parameters ---
num_fermions = 64 #! cost
dimensions = 2
bead_counts = [2]
tau_start = 5.25
tau_stop = 15.25 #! cost
tau_values = range(tau_start, tau_stop, length=48) #! cost

# --- High-Precision Settings ---
bigfloat_precision = 5000 #! cost
setprecision(ArbFloat, bigfloat_precision)

#Propagator choice: PA, TI, 4A is written as FA
propagator_choice = Symbol("PA")

# Now both thermo and hamiltonian estimators are calculated

# Use ArbFloat explicitly to prevent type promotion errors later
w = ArbFloat(1.0) 
balanced = false

run_id_T = "comp_thermo _N$num_fermions _D$dimensions _$propagator_choice"
run_id_H = "comp_ham _N$num_fermions _D$dimensions _$propagator_choice"