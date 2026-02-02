using ArbNumerics

# --- Simulation Parameters ---
num_fermions = 288
dimensions = 2
bead_counts = [200]
tau_start = 5.25
tau_stop = 15.25
tau_values = range(tau_start, tau_stop, length=16)

# --- High-Precision Settings ---
bigfloat_precision = 96000
setprecision(ArbFloat, bigfloat_precision)

#Propagator choice: PA, TI, 4A is written as FA
propagator_choice = Symbol("PA")

# Now both thermo and hamiltonian estimators are calculated

# Use ArbFloat explicitly to prevent type promotion errors later
w = ArbFloat(1) 
balanced = false

run_id_T = "comp_thermo _N$num_fermions _D$dimensions _$propagator_choice"
run_id_H = "comp_ham _N$num_fermions _D$dimensions _$propagator_choice"

#Preview plot
toplim = 20000
bottomlim = 5000