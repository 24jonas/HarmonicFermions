module HeatC

using CSV
using DataFrames

#include("parameters_ARB.jl")

function Get_heat_capacity()
    # --- Simulation Parameters ---
    num_fermions = 1024
    dimensions = 2
    propagator_choice = Symbol("PA")

    run_id_T = "comp_thermo _N$num_fermions _D$dimensions _$propagator_choice"
    run_id_H = "comp_ham _N$num_fermions _D$dimensions _$propagator_choice"

    data_file_T = "data_$(run_id_T).csv"
    data_file_H = "data_$(run_id_H).csv"

    # Read the data
    df_T = CSV.read(data_file_T, DataFrame)
    df_H = CSV.read(data_file_H, DataFrame)

    #--- Calculate Heat Capacity ---
    # Calculate the midpoints for the tau array (reduces size by exactly 1)
    n_tau = length(df_T.tau)
    midpoint_taus = [(df_T.tau[i] + df_T.tau[i+1]) / 2 for i in 1:(n_tau - 1)]
    Temp = midpoint_taus .^ (-1) # Temperature is the inverse of tau

    # Create new DataFrames using the midpoint tau column
    results_df_T = DataFrame(tau = Temp)
    results_df_H = DataFrame(tau = Temp)

    # Helper function to compute the discrete derivative 
    # C = -tau^2 * (dE / dtau)
    function compute_heat_capacity(tau_vals, E_vals)
        n = length(tau_vals)
        C_v = zeros(n - 1) # Array is exactly 1 element shorter
        
        for i in 1:(n - 1)
            # Pairwise difference (dE / dtau)
            dE_dtau = (E_vals[i+1] - E_vals[i]) / (tau_vals[i+1] - tau_vals[i])
            
            # The corresponding tau is the midpoint
            tau_mid = (tau_vals[i+1] + tau_vals[i]) / 2
            
            C_v[i] = -(tau_mid^2) * dE_dtau
        end
        
        return C_v
    end

    # Apply the computation to all N_bead columns in the Thermodynamics dataframe
    for col in names(df_T)
        if col != "tau"
            results_df_T[!, col] = compute_heat_capacity(df_T.tau, df_T[!, col])/num_fermions
        end
    end

    # Apply the computation to all N_bead columns in the Hamiltonian dataframe
    for col in names(df_H)
        if col != "tau"
            results_df_H[!, col] = compute_heat_capacity(df_H.tau, df_H[!, col])/num_fermions
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