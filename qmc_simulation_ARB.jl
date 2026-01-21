module QMC

using Polynomials
using DataFrames
using CSV
using JSON

include("parameters.jl")
include("PropagatorsModule.jl")
using .PropagatorsModule

include("EnergySol.jl")
using .EnergySol

using ArbNumerics


function run_and_plot()
    
    # --- Save Parameters ---

    # Explicitly set precision again to be safe
    setprecision(ArbFloat, bigfloat_precision)
    print("... set precision to $bigfloat_precision bits \n")

    p_funcs = propagators[propagator_choice]

    results_df = DataFrame(tau = collect(tau_values))

    factor_calc = get_factor(mode)
    energy_calc = get_energy_calc(balanced)

    print("... plot set up \n")

    # --- Main Loop ---
    for N in bead_counts
        energies = Vector{Float64}(undef, length(tau_values)) 

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
            gamma_val_s = (sqrt(zeta_1_s^2-1))/p_funcs.k1(epsilon)
            
            u_s = (zeta_1_s >= 1) ? acosh(zeta_1_s) : ArbFloat(0)
            b_s = exp(-N * u_s)
            
            # --- Energy Calculation ---
            # IMPORTANT: We do not use threading or recursion here.
            energy1 = harmEnergy(1, dimensions, b)
            energy1star = harmEnergy(1, dimensions, b_s)
            
            # Pass 'nothing' for cache compatibility
            energystar = energy_calc(num_fermions, dimensions, b_s, nothing)

            factor_regular, factor_star = factor_calc(lambda_val, gamma_val, w, lambda_val_s, gamma_val_s)

            energy = factor_regular*energy1 + factor_star*(energystar - energy1star) 
            
            energies[i] = energy
            
            # MANUAL GC:
            # 230,000 bits is ~30KB per number.
            # A single iteration creates megabytes of temporary garbage.
            # We force a cleanup every few steps to keep RAM usage stable.
            if i % 2 == 0
                GC.gc()
            end
        end

        results_df[!, "N_$(N)"] = energies

        # Save partial data
        CSV.write("partial_data_$(run_id).csv", results_df)
    end

    # --- Final Save ---
    output_filename = "data_$(run_id).csv"
    CSV.write(output_filename, results_df)
    display(plt)
end

end # module

QMC.run_and_plot()