using ArbNumerics

# --- Simulation Parameters ---
num_fermions = 1024 #! cost
dimensions = 2
bead_counts = [200] #! cost
tau_start = 5.00
tau_stop = 15.00 #! cost
tau_values = range(tau_start, tau_stop, length=16) #! resolution (proportional cost)

# tau_beads = [round(Int, tau * 4) for tau in tau_values] # for variable N program

# --- High-Precision Settings ---
bigfloat_precision = 256 #! cost, determined by n, tau and to lesser extent N
setprecision(ArbFloat, bigfloat_precision)

#Propagator choice: PA, TI, 4A is written as FA
propagator_choice = Symbol("PA")

# Now both thermo and hamiltonian estimators are calculated

# Use ArbFloat explicitly to prevent type promotion errors later
w = ArbFloat(0.5) 
balanced = false

run_id_T = "comp_thermo _N$num_fermions _D$dimensions _$propagator_choice"
run_id_H = "comp_ham _N$num_fermions _D$dimensions _$propagator_choice"