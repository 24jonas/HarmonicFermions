module HeatC

using CSV
using DataFrames

#include("parameters_ARB.jl")

function Get_heat_capacity()
    # --- Simulation Parameters ---
    num_fermions = 8
    dimensions = 1
    propagator_choice = Symbol("PA")

    run_id_T = "comp_thermo _N$num_fermions _D$dimensions _$propagator_choice"
    run_id_H = "comp_ham _N$num_fermions _D$dimensions _$propagator_choice"

    data_file_T = "data_$(run_id_T).csv"
    data_file_H = "data_$(run_id_H).csv"

    # Read the data
    df_T = CSV.read(data_file_T, DataFrame)
    df_H = CSV.read(data_file_H, DataFrame)

    #--- Calculate Heat Capacity ---
    # Create new DataFrames for the results, starting with the tau column
    results_df_T = DataFrame(tau = df_T.tau)
    results_df_H = DataFrame(tau = df_H.tau)

    # Helper function to compute the discrete derivative 
    # C = -tau^2 * (dE / dtau)
    function compute_heat_capacity(tau_vals, E_vals)
        n = length(tau_vals)
        C_v = zeros(n)
        for i in 1:n
            if i == 1
                # Forward difference for the first element
                dE_dtau = (E_vals[2] - E_vals[1]) / (tau_vals[2] - tau_vals[1])
            elseif i == n
                # Backward difference for the last element
                dE_dtau = (E_vals[n] - E_vals[n-1]) / (tau_vals[n] - tau_vals[n-1])
            else
                # Central difference for interior points
                dE_dtau = (E_vals[i+1] - E_vals[i-1]) / (tau_vals[i+1] - tau_vals[i-1])
            end
            C_v[i] = -(tau_vals[i]^2) * dE_dtau
        end
        return C_v
    end

    # Apply the computation to all N_bead columns in the Thermodynamics dataframe
    for col in names(df_T)
        if col != "tau"
            results_df_T[!, col] = compute_heat_capacity(df_T.tau, df_T[!, col])
        end
    end

    # Apply the computation to all N_bead columns in the Hamiltonian dataframe
    for col in names(df_H)
        if col != "tau"
            results_df_H[!, col] = compute_heat_capacity(df_H.tau, df_H[!, col])
        end
    end

    # Save
    # Saving with a _Cv suffix to avoid overwriting your raw energy data files
    output_filename_T = "data_Cv_$(run_id_T).csv"
    output_filename_H = "data_Cv_$(run_id_H).csv"
    
    CSV.write(output_filename_T, results_df_T)
    CSV.write(output_filename_H, results_df_H)
    
    println("Heat capacity computations saved successfully!")
end

end # module

HeatC.Get_heat_capacity()