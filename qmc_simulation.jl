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
        "bigfloat_precision" => bigfloat_precision
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

    # --- Main Loop ---
    for N in bead_counts
        # CHANGE 1: Pre-allocate the results array. `undef` is fine since we'll fill every spot.
        energies = Vector{Float64}(undef, length(tau_values)) 
        
        # CHANGE 2: Add the `@threads` macro to parallelize this loop.
        # We loop over indices (1, 2, 3...) instead of values to ensure thread safety.
        @threads for i in eachindex(tau_values)
            tau = tau_values[i]
    
            # --- The calculation logic remains identical ---
            epsilon = BigFloat(tau) / BigFloat(N)
            zeta_1 = p_funcs.zeta_1(epsilon)
            lambda_val = p_funcs.lambda(epsilon)
            gamma_val = p_funcs.gamma(epsilon)
            u = (zeta_1 >= 1) ? acosh(zeta_1) : BigFloat(0)
            b = exp(-N * u) #exp(-tau) #exp(-N * u)

            # --- Effective case ---
            zeta_1_s = p_funcs.zeta_1(w*epsilon)
            lambda_val_s = p_funcs.lambda(epsilon)
            gamma_val_s = p_funcs.gamma(epsilon)
            u_s = (zeta_1_s >= 1) ? acosh(zeta_1_s) : BigFloat(0)
            b_s = exp(-N * u_s) #exp(-tau) #exp(-N * u_s)
            
            # --- NEW: Call the direct-space recursive function ---
            cache = Dict{Int, Tuple{BigFloat, BigFloat}}()
            energy1 = harmEnergy(1, dimensions, b, cache)

            cache2 = Dict{Int, Tuple{BigFloat, BigFloat}}()
            cache3 = Dict{Int, Tuple{BigFloat, BigFloat}}()
            energy1star = harmEnergy(1, dimensions, b_s, cache2)

            if balanced
                energystar = harmEnergy(floor(Int,num_fermions/2), dimensions, b_s, cache2) + harmEnergy(ceil(Int,num_fermions/2), dimensions, b_s, cache3)
            else
                energystar = harmEnergy(num_fermions, dimensions, b_s, cache3)
            end


            if mode == "thermo"
                factor_regular = (lambda_val / gamma_val) #(1.0 / sqrt(1 + (epsilon^2) /4)) #PA
                factor_star = (w / sqrt(1 + (w^2) * (epsilon^2) / 4)) #PA
            else
                factor_regular = (gamma_val^2 + 1 )/(2*lambda_val)*(lambda_val / gamma_val) #0.5*(sqrt(1+(epsilon^2)/4) + 1.0 / sqrt(1 + (epsilon^2) /4)) #PA
                factor_star = (w*0.5)* (sqrt(1 + (epsilon^2)*(w^2)/4)+ 1.0 / sqrt(1 + (w^2) * (epsilon^2) / 4)) #PA
            end

            #cache0 = Dict{Int, Tuple{BigFloat, BigFloat}}()
            #energy_n2c = harmEnergy(ceil(Int,num_fermions/2), dimensions, b, cache0)
            #alt bal factor_regular*(energy1*0+energy_n2c) + factor_star*(energy1star)

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