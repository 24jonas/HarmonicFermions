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
        "tau_range" => Dict("start" => tau_start, "stop" => tau_stop),
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
    p_funcs = propagators[propagator_choice]

    # --- Plotting Setup ---
    # Discrete derivative gives energy at the interval midpoints
    tau_midpoints = (tau_values[1:end-1] .+ tau_values[2:end]) ./ 2.0
    results_df = DataFrame(tau = tau_midpoints)

    plt = plot(
        title = "High-Precision Energy vs. Tau (Discrete Derivative)",
        xlabel = "τ (Imaginary Time)",
        ylabel = "Energy E",
        legend = :topright,
        ylim = (bottomlim, toplim),
        grid = true
    )

    factor_calc = get_factor(mode)

    # --- Main Loop ---
    for N in bead_counts
        println("Processing N = $N...")

        # 1. Pre-calculate b values (Thread-safe)
        b_values = Vector{BigFloat}(undef, length(tau_values))
        b_values_s = Vector{BigFloat}(undef, length(tau_values))

        @threads for i in eachindex(tau_values)
            tau = tau_values[i]
            epsilon = BigFloat(tau) / BigFloat(N)
            
            # --- Standard Propagator (w=1) ---
            zeta_1 = p_funcs.zeta_1(epsilon)
            u = (zeta_1 >= 1) ? acosh(zeta_1) : BigFloat(0)
            b_values[i] = exp(-N * u)

            # --- Effective Propagator (Interaction w) ---
            zeta_1_s = p_funcs.zeta_1(w*epsilon)
            u_s = (zeta_1_s >= 1) ? acosh(zeta_1_s) : BigFloat(0)
            b_values_s[i] = exp(-N * u_s)
        end

        # 2. Compute Energies via Discrete Derivative
        # IMPORTANT: 'energy_1_std' must be for 1 particle (COM) to match the formula logic.
        
        # A. Center of Mass Energy at w=1 (1 Particle)
        _, energy_1_std = EnergySol.calculate_thermo_energy_discrete(tau_values, b_values, 1, dimensions)
        
        # B. Total Energy at w (N Particles)
        _, energy_tot_star = EnergySol.calculate_thermo_energy_discrete(tau_values, b_values_s, num_fermions, dimensions)
        
        # C. Center of Mass Energy at w (1 Particle)
        _, energy_1_star = EnergySol.calculate_thermo_energy_discrete(tau_values, b_values_s, 1, dimensions)

        # 3. Apply the COM Correction Formula
        # Formula: E_final = E_COM(w=1) + [ E_Total(w) - E_COM(w) ]
        # This replaces the Relative Energy at w=1 with Relative Energy at w, 
        # while keeping the COM Energy at w=1.
        
        final_energies = Vector{Float64}(undef, length(energy_1_std))

        for k in eachindex(energy_1_std)
            
            # Note: Since calculate_thermo_energy_discrete returns the Physical Energy directly,
            # the implicit factors are 1.0. Applying analytic factors (like w^2) here would 
            # incorrectly scale the already-physical energies. We preserve the STRUCTURE of 
            # the user's formula which performs the physics separation.
            
            term_com = energy_1_std[k]               # E_COM(w=1)
            term_rel = energy_tot_star[k] - energy_1_star[k] # E_Rel(w)
            
            final_energies[k] = term_com + term_rel
        end

        results_df[!, "N_$(N)"] = final_energies
        plot!(plt, tau_midpoints, final_energies, label="N = $N")
    end

    # --- Save and Display ---
    output_filename = "data_$(run_id).csv"
    CSV.write(output_filename, results_df)

    display(plt)
end

end # module QMC

# --- To Run the Code ---
QMC.run_and_plot()