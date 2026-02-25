using ArbNumerics

# --- Simulation Parameters ---
num_fermions = 8 #! cost
dimensions = 1
bead_counts = [2,4,6,8,10,12,14,16,18,20,22,24,26,28,30,32,34,36,38,40,42,44,46,48,50,52,54,56,58,60,62,64,66,68,70,72,74,76,78,80,82,84,86,88,90,92,94,96,98,100] #! cost
tau_start = 0.25
tau_stop = 6.25 #! cost
tau_values = range(tau_start, tau_stop, length=32) #! resolution (proportional cost)

# --- High-Precision Settings ---
bigfloat_precision = 5000 #! cost, determined by n, tau and to lesser extent N
setprecision(ArbFloat, bigfloat_precision)

#Propagator choice: PA, TI, 4A is written as FA
propagator_choice = Symbol("PA")

# Now both thermo and hamiltonian estimators are calculated

# Use ArbFloat explicitly to prevent type promotion errors later
w = ArbFloat(1.0) 
balanced = false

run_id_T = "comp_thermo _N$num_fermions _D$dimensions _$propagator_choice"
run_id_H = "comp_ham _N$num_fermions _D$dimensions _$propagator_choice"