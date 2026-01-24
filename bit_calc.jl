module PrecisionCheck

using ArbNumerics
using Printf

# --- LOAD YOUR MODULES ---
include("PropagatorsModule.jl")
using .PropagatorsModule

# --- SETUP ---
num_fermions = 500      # FIXED: The number of particles (Recursion Depth)
dimensions = 2
beta = 15.25            # Max imaginary time (Worst case)
bead_counts = [200] # We check scaling across these bead counts

# Toggle this to see the difference
propagator_choice = Symbol("PA") 
# propagator_choice = Symbol("FA") 

# Low precision is fine for the *estimate* (we just need the exponent)
setprecision(BigFloat, 2048) 

function check_scaling(num_fermions, D, beta, bead_list, prop_choice)
    println("--- PRECISION SCALING CHECK ---")
    println("Particles: $num_fermions | Beta: $beta | Propagator: $prop_choice")
    println("-------------------------------------------------------------")
    @printf "%-10s | %-15s | %-15s | %-10s\n" "Beads" "Log(Noise)" "Log(Signal)" "Bits Needed"
    println("-------------------------------------------------------------")

    # 1. CALCULATE FIXED FERMIONIC SIGNAL (Physical Constant)
    # This does not depend on beads, only on particle count/beta.
    E_Fermion = 0.0
    count = 0
    k = 1
    while count < num_fermions
        degeneracy = k 
        take = min(degeneracy, num_fermions - count)
        E_Fermion += take * k 
        count += take
        k += 1
    end
    log_Z_F = -beta * E_Fermion # The Signal

    # 2. LOOP OVER BEAD COUNTS
    for P in bead_list
        # P is the number of beads (Time Slices)
        
        # --- A. Get b from Propagator ---
        p_funcs = propagators[prop_choice]
        epsilon = ArbFloat(beta) / ArbFloat(P) # Epsilon depends on BEADS
        
        zeta_1 = p_funcs.zeta_1(epsilon)
        u = (zeta_1 >= 1) ? acosh(zeta_1) : ArbFloat(0)
        b = exp(-P * u)
        
        b_val = BigFloat(b)

        # --- B. Calculate Z_Boson (Noise) ---
        # The sums run up to num_fermions (Particles)
        z_vals = Vector{BigFloat}(undef, num_fermions)
        b_pow = b_val
        
        for i in 1:num_fermions
            num = b_pow^(D/2)
            den = (1 - b_pow)^D
            z_vals[i] = abs(num / den)
            if i < num_fermions; b_pow *= b_val; end
        end
        
        Z = Vector{BigFloat}(undef, num_fermions + 1)
        Z[1] = 1.0
        
        for k in 1:num_fermions
            sum_Z = BigFloat(0.0)
            for i in 1:k
                sum_Z += z_vals[i] * Z[k - i + 1]
            end
            Z[k+1] = sum_Z / k
        end
        
        log_Z_B = log(Z[num_fermions+1])
        
        # --- C. Compare ---
        bits_needed = (log_Z_B - log_Z_F) / log(2)
        
        @printf "%-10d | %-15.2f | %-15.2f | %-10.0f\n" P log_Z_B log_Z_F bits_needed
    end
    println("-------------------------------------------------------------")
    println("Note: Add ~64-128 bits buffer to 'Bits Needed' for safety.")
end

check_scaling(num_fermions, dimensions, beta, bead_counts, propagator_choice)

end # module