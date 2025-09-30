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

# Ensure you have a plotting backend, e.g., Pkg.add("GR")
gr() 

#=============================================================================
 Part 1: Exact Ground State Energy
=============================================================================#


#=============================================================================
 Part 2: High-Precision Coupled Recursion in Log-Space
=============================================================================#

function z_and_derivative_high_precision(i::Int, d::Int, b::BigFloat)
    if !(0 < b < 1)
        return (zero(BigFloat), zero(BigFloat))
    end
    
    b_i = b^i
    one_minus_b_i = 1 - b_i
    
    # Direct formula for z_i = b^(i*d/2) / (1 - b^i)^d
    z = b^(BigFloat(i * d) / 2) / (one_minus_b_i)^d
    
    # We still need the logarithmic derivative to find z_i'
    # (log(z))' = (i*d)/(2*b) + (d*i*b^(i-1))/(1-b^i)
    log_deriv = (BigFloat(i * d) / (2 * b)) + (d * i * b^(i - 1)) / one_minus_b_i
    
    # Direct formula for the derivative: z_i' = z_i * (log(z_i))'
    z_prime = z * log_deriv
    
    return z, z_prime
end

function Z_and_derivative_recursive!(n::Int, d::Int, b::BigFloat, cache::Dict)
    # The cache now stores (Z_n, Z_n') instead of their logs.
    if haskey(cache, n)
        return cache[n]
    end
    # Base case: Z_0 = 1, Z_0' = 0
    if n == 0
        return (one(BigFloat), zero(BigFloat))
    end

    total_Z = zero(BigFloat)
    total_Z_prime = zero(BigFloat)

    for i in 1:n
        # --- Get z_i and its derivative z_i' ---
        z_i, z_i_prime = z_and_derivative_high_precision(i, d, b)

        # --- Recursive call for the remaining n-i particles ---
        Z_ni, Z_ni_prime = Z_and_derivative_recursive!(n - i, d, b, cache)
        
        sign = iseven(i - 1) ? 1 : -1

        # --- Sum the terms directly ---
        # For Z_n
        total_Z += sign * z_i * Z_ni
        
        # For Z_n', using the product rule: (f*g)' = f'*g + f*g'
        total_Z_prime += sign * (z_i_prime * Z_ni + z_i * Z_ni_prime)
    end

    # Final result is (1/n) * sum
    final_Z = total_Z / n
    final_Z_prime = total_Z_prime / n
    
    cache[n] = (final_Z, final_Z_prime)
    return final_Z, final_Z_prime
end


#=============================================================================
 Part 3: Main Plotting Script
=============================================================================#

function run_and_plot()
    # Load all parameters from the external file
    include("parameters.jl")

    params_to_save = Dict(
        "mode" => mode,
        "run_id" => run_id,
        "num_fermions" => num_fermions,
        "dimensions" => dimensions,
        "propagator_choice" => propagator_choice,
        "bead_counts" => bead_counts,
        "tau_range" => Dict(
            "start" => tau_start,
            "stop" => tau_stop,
            "length" => tau_length
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
    propagators = (
        # Primitive Action (2nd Order, Error ~ ε³) - what you had before
        PA = (
            zeta_1 = e -> 1 + e^2/2,
            lambda = e -> BigFloat(1.0),
            gamma = e -> sqrt(1 + e^2/4)
        ),

        # Takahashi-Imada Propagator (4th Order, Error ~ ε⁵)
        TI = (
            zeta_1 = e -> begin
                omega_eff_sq = 1 + e^2 / 12
                omega_eff = sqrt(omega_eff_sq)
                return cosh(omega_eff * e)
            end,
            lambda = e -> BigFloat(1.0), # For the harmonic oscillator, lambda is simple
            gamma = e -> begin
                omega_eff_sq = 1 + e^2 / 12
                omega_eff = sqrt(omega_eff_sq)
                # Use sinc function for numerical stability near e=0
                return sqrt(sinh(omega_eff * e) / (omega_eff * e))
            end
        ),

        # A Fourth-Order Propagator (4A, Error ~ ε⁵) - based on work by S. A. Chin
        # This provides an alternative high-accuracy method.
        FA = (
            zeta_1 = e -> (1 + e^2/6) * cosh(e) - e^2/6,
            lambda = e -> 1 - e^2/12,
            gamma = e -> sqrt( ((1+e^2/6)*sinh(e) - (e^3/6)*cosh(e)) / e )
        )
    )
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
            # The cache is now a Dict for (BigFloat, BigFloat)
            cache = Dict{Int, Tuple{BigFloat, BigFloat}}()
            
            # --- The calculation logic remains identical ---
            epsilon = BigFloat(tau) / BigFloat(N)
            zeta_1 = p_funcs.zeta_1(epsilon)
            lambda_val = p_funcs.lambda(epsilon)
            gamma_val = p_funcs.gamma(epsilon)
            u = (zeta_1 >= 1) ? acosh(zeta_1) : BigFloat(0)
            b = exp(-N * u) #exp(-tau) #exp(-N * u)
            
            # --- NEW: Call the direct-space recursive function ---
            Z_val, Z_prime = Z_and_derivative_recursive!(
                1, dimensions, b, cache
            )

            # --- NEW: Update the energy calculation ---
            # The logarithmic derivative is Z'/Z.
            energy1 = if isnan(Z_prime) || Z_val <= 0 # Avoid division by zero or log of negative
                NaN
            else
                log_Z_prime = Z_prime / Z_val # This is the log derivative
                # The rest of the formula is the same
                Float64((b * log_Z_prime))
            end

            w = BigFloat(3/2)
            cache2 = Dict()
            cache3 = Dict()

            # --- The calculation logic remains identical ---
            epsilon = BigFloat(tau) / BigFloat(N)
            zeta_1 = p_funcs.zeta_1(w*epsilon)
            lambda_val = p_funcs.lambda(epsilon)
            gamma_val = p_funcs.gamma(epsilon)
            u = (zeta_1 >= 1) ? acosh(zeta_1) : BigFloat(0)
            b = exp(-N * u) #exp(-tau) #exp(-N * u)

            # --- NEW: Call the direct-space recursive function ---
            Z1_val, Z1_prime = Z_and_derivative_recursive!(
                1, dimensions, b, cache2
            )

            # --- NEW: Call the direct-space recursive function ---
            Z_val_f, Z_prime_f = Z_and_derivative_recursive!(
                floor(Int,num_fermions/2), dimensions, b, cache2
            )

            # --- NEW: Call the direct-space recursive function ---
            Z_val_c, Z_prime_c = Z_and_derivative_recursive!(
                ceil(Int,num_fermions/2), dimensions, b, cache3
            )


            energy1star = if isnan(Z1_prime) || Z1_val <= 0 # Avoid division by zero or log of negative
                NaN
            else
                log_Z_prime = Z1_prime / Z1_val # This is the log derivative
                # The rest of the formula is the same
                Float64((b * log_Z_prime))
            end

            energystar_c = if isnan(Z_prime_c) || Z_val_c <= 0 # Avoid division by zero or log of negative
                NaN
            else
                log_Z_prime = Z_prime_c / Z_val_c # This is the log derivative
                # The rest of the formula is the same
                Float64((b * log_Z_prime))
            end

            energystar_f = if isnan(Z_prime_f) || Z_val_f <= 0 # Avoid division by zero or log of negative
                NaN
            else
                log_Z_prime = Z_prime_f / Z_val_f # This is the log derivative
                # The rest of the formula is the same
                Float64((b * log_Z_prime))
            end

            if mode == "thermo"
                energy_c = (1.0 / sqrt(1 + (epsilon^2) /4))*energy1 + (w / sqrt(1 + (w^2) * (epsilon^2) / 4))*(energystar_c - energy1star)
            else
                energy_c = 0.5*(sqrt(1+(epsilon^2)/4) + 1.0 / sqrt(1 + (epsilon^2) /4))*energy1 + (w*0.5)* (sqrt(1 + (epsilon^2)*(w^2)/4)+ 1.0 / sqrt(1 + (w^2) * (epsilon^2) / 4))*(energystar_c - energy1star)
            end

            if mode == "thermo"
                energy_f = (1.0 / sqrt(1 + (epsilon^2) /4))*energy1 + (w / sqrt(1 + (w^2) * (epsilon^2) / 4))*(energystar_f - energy1star)
            else
                energy_f = 0.5*(sqrt(1+(epsilon^2)/4) + 1.0 / sqrt(1 + (epsilon^2) /4))*energy1 + (w*0.5)* (sqrt(1 + (epsilon^2)*(w^2)/4)+ 1.0 / sqrt(1 + (w^2) * (epsilon^2) / 4))*(energystar_f - energy1star)
            end
            
            energies[i] = energy_c + energy_f
        end

        results_df[!, "N_$(N)"] = energies

        plot!(plt, tau_values, energies, label="N = $N beads")
    end

    # --- Add Exact Energy Line ---
    #exact_energy = get_ground_state_energy(num_fermions, dimensions)
	#exact_energy = num_fermions^2 / 2.0
    #hline!(plt, [exact_energy], linestyle=:dash, color=:black, label="Exact Energy = $(round(exact_energy, digits=2))")

    # --- Display the Plot ---
    output_filename = "data_$(run_id).csv"
    CSV.write(output_filename, results_df)

    display(plt)
end

end # module QMC

# --- To Run the Code ---
QMC.run_and_plot()