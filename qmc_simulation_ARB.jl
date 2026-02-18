module QMC

using Polynomials
using DataFrames
using CSV
using JSON

include("parameters_ARB.jl")
include("PropagatorsModule_ARB.jl")
using .PropagatorsModule

include("EnergySol_ARB.jl")
using .EnergySol

using Base.Threads

using ArbNumerics


function run_and_plot()

    # ========================================================
    # MEMORY USAGE DIAGNOSTIC (OPTIMIZED IMPLEMENTATION)
    # ========================================================

    # ========================================================
    # MEMORY USAGE DIAGNOSTIC (CORRECTED FOR C-HEAP)
    # ========================================================

    # Explicitly set precision again to be safe
    setprecision(ArbFloat, bigfloat_precision)
    print("... set precision to $bigfloat_precision bits \n")

    println("------------------------------------------------")
    println("Memory Diagnostic (N=$num_fermions, Threads=$(Threads.nthreads()), Precision=$(precision(ArbFloat)) bits):")

    # 1. Calculate Theoretical Size of One ArbFloat
    # ArbFloat stores: 
    #   - A mid-point (mantissa) -> requires 'precision' bits
    #   - A radius (error bound) -> relatively small, fixed size mag_t
    #   - Struct overhead -> ~32-64 bytes
    
    bits = precision(ArbFloat)
    bytes_per_number = (bits / 8) + 64 # +64 bytes for struct/radius overhead
    
    # 2. Active Working Memory (Per Thread)
    # The new 'harmEnergy_fast' uses 4 vectors simultaneously:
    # z_vals, z_prime_vals, Z, Z_prime
    # Each vector has 'num_fermions' elements.
    
    total_elements_per_thread = num_fermions * 4
    mem_per_thread_bytes = total_elements_per_thread * bytes_per_number
    mem_total_active_bytes = mem_per_thread_bytes * Threads.nthreads()

    # 3. Permanent Arrays (Output)
    # This is for 'energies' (Float64), which is trivial. 
    # If using 'calculate_thermo_energy_discrete', you might store b_values or Z_values.
    mem_permanent_bytes = 0.0
    
    # Convert to MB
    mb_active = mem_total_active_bytes / 1024^2
    mb_per_thread = mem_per_thread_bytes / 1024^2

    println("  Est. Size per Number: ~$(round(bytes_per_number/1024, digits=4)) KB")
    println("  Per-Thread Cache:     ~$(round(mb_per_thread, digits=4)) MB")
    println("  TOTAL Active Memory:  ~$(round(mb_active, digits=4)) MB")
    println("------------------------------------------------")
    
    # --- Save Parameters ---

    p_funcs = propagators[propagator_choice]

    results_df_T = DataFrame(tau = collect(tau_values))
    results_df_H = DataFrame(tau = collect(tau_values))

    factor_calc_T = get_factor("thermo")
    factor_calc_H = get_factor("ham")
    energy_calc = get_energy_calc(balanced)

    # --- Main Loop ---
    for N in bead_counts
        energies_T = Vector{Float64}(undef, length(tau_values))
        energies_H = Vector{Float64}(undef, length(tau_values)) 

        print("... working on bead $N \n")
        
        # SINGLE THREADED EXECUTION to avoid stack/memory collisions with huge types
        for i in eachindex(tau_values)
            tau = tau_values[i]
            
            # Use ArbFloat explicitly for all conversions
            epsilon = ArbFloat(tau) / ArbFloat(N)
            
            zeta_1 = p_funcs.zeta_1(epsilon)
            lambda_val = p_funcs.lambda(epsilon)
            gamma_val = p_funcs.gamma(epsilon) 

            u = (zeta_1 >= 1) ? acosh(zeta_1) : ArbFloat(0)
            b = exp(-N * u)

            # --- Effective case ---
            zeta_1_s = p_funcs.zeta_1(w*epsilon)
            lambda_val_s = p_funcs.lambda(w*epsilon) 
            gamma_val_s = (sqrt(zeta_1_s^2-1))/p_funcs.k1(w*epsilon)
            
            u_s = (zeta_1_s >= 1) ? acosh(zeta_1_s) : ArbFloat(0)
            b_s = exp(-N * u_s)
            
            # --- Energy Calculation ---
            # IMPORTANT: We do not use threading or recursion here.
            energy1 = harmEnergy(1, dimensions, b)
            energy1star = harmEnergy(1, dimensions, b_s)
            
            # Pass 'nothing' for cache compatibility
            energystar = energy_calc(num_fermions, dimensions, b_s, nothing)

            factor_regular_T, factor_star_T = factor_calc_T(lambda_val, gamma_val, w, lambda_val_s, gamma_val_s)
            factor_regular_H, factor_star_H = factor_calc_H(lambda_val, gamma_val, w, lambda_val_s, gamma_val_s)

            energy_T = factor_regular_T*energy1 + factor_star_T*(energystar - energy1star) 
            energy_H = factor_regular_H*energy1 + factor_star_H*(energystar - energy1star) 

            energies_T[i] = energy_T
            energies_H[i] = energy_H
            
            # MANUAL GC:
            # 230,000 bits is ~30KB per number.
            # A single iteration creates megabytes of temporary garbage.
            # We force a cleanup every few steps to keep RAM usage stable.
            if i % 2 == 0
                GC.gc()
            end
        end

        results_df_T[!, "N_$(N)"] = energies_T
        results_df_H[!, "N_$(N)"] = energies_H

        # Save partial data
        CSV.write("partial_data_$(run_id_T).csv", results_df_T)
        CSV.write("partial_data_$(run_id_H).csv", results_df_H)
    end

    # --- Final Save ---
    output_filename_T = "data_$(run_id_T).csv"
    output_filename_H = "data_$(run_id_H).csv"
    CSV.write(output_filename_T, results_df_T)
    CSV.write(output_filename_H, results_df_H)
    print(results_df_H)
end

end # module





QMC.run_and_plot()