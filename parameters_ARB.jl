using ArbNumerics

# --- Simulation Parameters ---
num_fermions = 10 #! cost
dimensions = 1
bead_counts = [2, 20, 200] #! cost
tau_start = 0.005
tau_stop = 4.0 #! cost
tau_values = range(tau_start, tau_stop, length=80) #! resolution (proportional cost)

tau_beads = [round(Int, tau * 4) for tau in tau_values] # for variable N program

# --- High-Precision Settings ---
bigfloat_precision = 30000 #! cost, determined by n, tau and to lesser extent N
setprecision(ArbFloat, bigfloat_precision)

#Propagator choice: PA, TI, 4A is written as FA
propagator_choice = Symbol("PA")

# Now both thermo and hamiltonian estimators are calculated

# Use ArbFloat explicitly to prevent type promotion errors later
w = ArbFloat(1.0) 
balanced = false

run_id_T = "comp_thermo _N$num_fermions _D$dimensions _$propagator_choice"
run_id_H = "comp_ham _N$num_fermions _D$dimensions _$propagator_choice"