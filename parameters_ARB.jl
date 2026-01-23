using ArbNumerics

# --- Simulation Parameters ---
num_fermions = 500
dimensions = 2
bead_counts = [2, 4, 8, 16]
tau_start = 0.25
tau_stop = 15.25
tau_values = range(tau_start, tau_stop, length=45)

# --- High-Precision Settings ---
# STRICT REQUIREMENT: 230,000 bits (~69,000 decimal digits)
bigfloat_precision = 230000 
setprecision(ArbFloat, bigfloat_precision)

#Propagator choice: PA, TI, 4A is written as FA
propagator_choice = Symbol("PA")

mode = "thermo" # thermo, ham

# Use ArbFloat explicitly to prevent type promotion errors later
w = ArbFloat(0.5) 
balanced = false

run_id = "comp_$mode _N$num_fermions _D$dimensions _$propagator_choice"

#Preview plot
toplim = 20000
bottomlim = 5000