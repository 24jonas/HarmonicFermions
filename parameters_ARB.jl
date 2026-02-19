using ArbNumerics

# --- Simulation Parameters ---
num_fermions = 96 #! cost
dimensions = 2
bead_counts = [96] #! cost
tau_start = 0.25
tau_stop = 12.25 #! cost
tau_values = range(tau_start, tau_stop, length=48) #! resolution (proportional cost)

# --- High-Precision Settings ---
bigfloat_precision = 6000 #! cost, determined by n, tau and to lesser extent N
setprecision(ArbFloat, bigfloat_precision)

#Propagator choice: PA, TI, 4A is written as FA
propagator_choice = Symbol("PA")

# Now both thermo and hamiltonian estimators are calculated

# Use ArbFloat explicitly to prevent type promotion errors later
w = ArbFloat(1.0) 
balanced = false

run_id_T = "comp_thermo _N$num_fermions _D$dimensions _$propagator_choice"
run_id_H = "comp_ham _N$num_fermions _D$dimensions _$propagator_choice"