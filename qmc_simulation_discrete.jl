# --- Add required packages ---
# In the Julia REPL, run:
# import Pkg
# Pkg.add(["Polynomials", "Plots"])

module QMC

using Polynomials
using Plots
using Base.Threads

using DataFrames
using CSV
using JSON

include("parameters.jl")

include("PropagatorsModule.jl")
using .PropagatorsModule

include("EnergySol.jl")
using .EnergySol

# Ensure you have a plotting backend, e.g., Pkg.add("GR")
gr() 


function run_and_plot()
    # Load all parameters from the external file
    
    params_to_save = Dict(
        "mode" => mode,
        "run_id" => run_id,
        "num_fermions" => num_fermions,
        "dimensions" => dimensions,
        "propagator_choice" => propagator_choice,
        "bead_counts" => bead_counts,
        "tau_range" => Dict(
            "start" => tau_start,
            "stop" => tau_stop
        ),
        "bigfloat_precision" => bigfloat_precision,
        "balanced" => balanced,
        "interaction w" => w
    )

    param_filename = "parameters_$(run_id).json"
    open(param_filename, "w") do f
        JSON.print(f, params_to_save, 4)
    end

	setprecision(BigFloat, bigfloat_precision)

    # --- Propagator Definition ---
    # Using a named tuple for clean access
    p_funcs = propagators[propagator_choice]

    # --- Plotting Setup ---
    plt = plot(
        title = "High-Precision Energy vs. Tau for n=$num_fermions, d=$dimensions",
        xlabel = "τ (Imaginary Time)",
        ylabel = "Energy E",
        legend = :topright,
        ylim = (bottomlim, toplim),
        grid = true
    )

    results_df = DataFrame(tau = collect(tau_values))

    factor_calc = get_factor(mode)
    energy_calc = get_energy_calc(balanced)

    # --- Main Loop ---
    for N in bead_counts
        # CHANGE 1: Pre-allocate the results array. `undef` is fine since we'll fill every spot.
        energies = Vector{Float64}(undef, length(tau_values)) 
        
        # CHANGE 2: Add the `@threads` macro to parallelize this loop.
        # We loop over indices (1, 2, 3...) instead of values to ensure thread safety.
        @threads for i in eachindex(tau_values)
            tau = tau_values[i]
    
            
            u = (zeta_1 >= 1) ? acosh(zeta_1) : BigFloat(0)
            # 1. Pre-calculate b values for all tau points
            b_values = Vector{BigFloat}(undef, length(tau_values))

            for i in eachindex(tau_values)
                tau = tau_values[i]
                epsilon = BigFloat(tau) / BigFloat(N)
                zeta_1 = p_funcs.zeta_1(epsilon)
                lambda_val = p_funcs.lambda(epsilon)
                gamma_val = p_funcs.gamma(epsilon) #p_funcs.gamma(epsilon) #goes to 1 for exact
                
                # Calculate propagator specifics (e.g., using PA)
                zeta_1 = p_funcs.zeta_1(epsilon)
                u = (zeta_1 >= 1) ? acosh(zeta_1) : BigFloat(0)
                b_values[i] = exp(-N * u)
            end



            # --- Effective case ---
            
            # 1. Pre-calculate b values for all tau points
            b_values_s = Vector{BigFloat}(undef, length(tau_values))

            for i in eachindex(tau_values)
                tau = tau_values[i]
                epsilon = BigFloat(tau) / BigFloat(N)
                zeta_1_s = p_funcs.zeta_1(w*epsilon)
                lambda_val_s = p_funcs.lambda(w*epsilon) 
                gamma_val_s = (sqrt(zeta_1_s^2-1))/p_funcs.k1(epsilon)#p_funcs.gamma(w*epsilon)#(sqrt(zeta_1_s^2-1))/p_funcs.k1(w*epsilon)
                u_s = (zeta_1_s >= 1) ? acosh(zeta_1_s) : BigFloat(0)
                
                b_values_s[i] = exp(-N * u_s)
            end

            # 2. Compute energies using the new discrete function
            # 2. Compute energies using the new discrete function
            tau_mid, energy = EnergySol.calculate_thermo_energy_discrete(tau_values, b_values, num_fermions, dimensions)
            tau_mid_s, energy1star = EnergySol.calculate_thermo_energy_discrete(tau_values, b_values_s, 1 dimensions)
            tau_mid_s, energystar = EnergySol.calculate_thermo_energy_discrete(tau_values, b_values_s, num_fermions, dimensions)

            factor_regular, factor_star = factor_calc(lambda_val, gamma_val, w, lambda_val_s, gamma_val_s)

            energy = factor_regular*energy1 + factor_star*(energystar -energy1star) 
            
            energies[i] = energy
        end

        results_df[!, "N_$(N)"] = energies

        plot!(plt, tau_values, energies, label="N = $N beads")
    end

    # --- Display the Plot ---
    output_filename = "data_$(run_id).csv"
    CSV.write(output_filename, results_df)

    display(plt)
end

end # module QMC

# --- To Run the Code ---
QMC.run_and_plot()