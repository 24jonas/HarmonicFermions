module PrecisionPredictorContinuum

using Printf

# --- PARAMETERS ---
N = 128             # Number of particles
beta = 15.25         # Inverse Temperature
dimensions = 2       # Dimensions (d=2)

function predict_bits_continuum(n_particles, d, beta_val)
    println("--- CONTINUUM PRECISION PREDICTOR ---")
    println("Particles: $n_particles | Dim: $d | Beta: $beta_val")
    println("---------------------------------------------")

    # 1. CALCULATE EXACT BOSON ENERGY (Ground State)
    # In the continuum limit (T -> 0), all Bosons sit in the lowest energy level.
    # 2D HO Ground State Energy per particle = d/2 = 1.0 (in units of hw)
    # Total E_B = N * E_0
    
    e_ground_single = d / 2.0
    E_Boson = n_particles * e_ground_single

    # 2. CALCULATE EXACT FERMION ENERGY (Stacked Shells)
    # We fill energy shells k=0, 1, 2... until we run out of particles.
    # Energy of shell k: E_k = k + d/2
    # Degeneracy of shell k: Binomial(k + d - 1, d - 1)
    
    E_Fermion = 0.0
    particles_left = n_particles
    k = 0 # Shell index (0 = ground state)

    while particles_left > 0
        # Calculate degeneracy for this shell
        if d == 1
            deg = 1
        elseif d == 2
            deg = k + 1
        else
            deg = 1
            for i in 1:(d-1)
                deg = deg * (k + i) // i
            end
        end

        # Fill the shell
        take = min(deg, particles_left)
        energy_of_shell = k + e_ground_single
        
        E_Fermion += take * energy_of_shell
        particles_left -= take
        k += 1
    end

    # 3. CALCULATE BITS NEEDED
    # Formula: Bits = Beta * (E_F - E_B) / ln(2)
    # This represents the entropy cost of the Pauli Exclusion Principle.
    
    energy_gap = E_Fermion - E_Boson
    bits_needed = (beta_val * energy_gap) / log(2)

    # --- OUTPUT ---
    @printf "Boson Energy (E_B):    %10.2f\n" E_Boson
    @printf "Fermion Energy (E_F):  %10.2f\n" E_Fermion
    @printf "Energy Gap (dE):       %10.2f\n" energy_gap
    println("---------------------------------------------")
    @printf "BITS NEEDED (Exact):   %10.0f\n" bits_needed
    @printf "Recommended (+128):    %10.0f\n" (bits_needed + 128)
end

predict_bits_continuum(N, dimensions, beta)

end #module