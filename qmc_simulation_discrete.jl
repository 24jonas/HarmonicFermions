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


# ========================================== Get peak memory
if Sys.iswindows()
    # --- WINDOWS IMPLEMENTATION ---
    struct PROCESS_MEMORY_COUNTERS
        cb::UInt32
        PageFaultCount::UInt32
        PeakWorkingSetSize::Csize_t  # Peak physical memory in Bytes
        WorkingSetSize::Csize_t
        QuotaPeakPagedPoolUsage::Csize_t
        QuotaPagedPoolUsage::Csize_t
        QuotaPeakNonPagedPoolUsage::Csize_t
        QuotaNonPagedPoolUsage::Csize_t
        PagefileUsage::Csize_t
        PeakPagefileUsage::Csize_t
    end

    function get_max_memory_usage_kb()
        # Get handle to current process
        hProcess = ccall(:GetCurrentProcess, stdcall, Ptr{Cvoid}, ())
        
        # Prepare struct
        pmc = Ref{PROCESS_MEMORY_COUNTERS}()
        # Set the size of the structure (cb)
        # Note: We must manually initialize the struct with 0s or just set cb after? 
        # Julia's Ref creates it. We just need to tell the API the size.
        cb = sizeof(PROCESS_MEMORY_COUNTERS)
        
        # Call Windows API (psapi.dll)
        ret = ccall((:GetProcessMemoryInfo, "psapi"), stdcall, Int32, 
                    (Ptr{Cvoid}, Ptr{PROCESS_MEMORY_COUNTERS}, UInt32), 
                    hProcess, pmc, cb)
        
        if ret == 0
            return -1 # Failed
        end
        
        # Windows reports in BYTES, convert to KB
        return Clong(pmc[].PeakWorkingSetSize / 1024)
    end

else
    # --- LINUX/MAC IMPLEMENTATION ---
    struct RUsage
        ru_utime_sec::Clong
        ru_utime_usec::Clong
        ru_stime_sec::Clong
        ru_stime_usec::Clong
        ru_maxrss::Clong            # Max resident set size
        ru_ixrss::Clong
        ru_idrss::Clong
        ru_isrss::Clong
        ru_minflt::Clong
        ru_majflt::Clong
        ru_nswap::Clong
        ru_inblock::Clong
        ru_oublock::Clong
        ru_msgsnd::Clong
        ru_msgrcv::Clong
        ru_nsignals::Clong
        ru_nvcsw::Clong
        ru_nivcsw::Clong
    end

    function get_max_memory_usage_kb()
        ru = Vector{RUsage}(undef, 1)
        ret = ccall(:getrusage, Cint, (Cint, Ptr{Cvoid}), 0, ru)
        if ret == 0
            val = ru[1].ru_maxrss
            # macOS reports in Bytes, Linux in KB
            if Sys.isapple()
                return Clong(val / 1024)
            else
                return val
            end
        else
            return -1
        end
    end
end
# ==========================================


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

        # ========================================================
        # MEMORY USAGE DIAGNOSTIC
        # ========================================================
        
        # 1. Permanent Arrays (Scaling with Time Steps)
        # We measure one vector (b_values) and estimate the total based on how many
        # similar vectors are active simultaneously (approx 5: b, b_s, Z_std, Z_star, Z_1star).
        mem_single_array = Base.summarysize(b_values)
        mem_permanent_arrays = mem_single_array * 5
        
        # 2. Active Cache (Scaling with Particles and Threads)
        # The cache is temporary, so we simulate one to measure its size.
        # It stores N BigFloats. We multiply by nthreads() because each thread has its own.
        dummy_cache = Dict{Int, BigFloat}()
        for k in 1:num_fermions
            dummy_cache[k] = BigFloat(1.0)
        end
        mem_single_cache = Base.summarysize(dummy_cache)
        mem_total_cache = mem_single_cache * Threads.nthreads()
        
        # Convert to MB
        mb_arrays = mem_permanent_arrays / 1024^2
        mb_cache = mem_total_cache / 1024^2
        
        println("------------------------------------------------")
        println("Memory Diagnostic (N=$N, Threads=$(Threads.nthreads())):")
        println("  Permanent Arrays (Timesteps x precision): ~$(round(mb_arrays, digits=6)) MB")
        println("  Active Cache (n x threads x precision): ~$(round(mb_cache, digits=6)) MB")

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

    # --- REPORT PEAK MEMORY ---
    max_mem_kb = get_max_memory_usage_kb()
    max_mem_mb = max_mem_kb / 1024.0
    max_mem_gb = max_mem_mb / 1024.0
    
    println("------------------------------------------------")
    println("Peak Memory Usage (MaxRSS):")
    println("  $(max_mem_kb) KB")
    println("  $(round(max_mem_mb, digits=2)) MB")
    println("  $(round(max_mem_gb, digits=3)) GB")
    println("------------------------------------------------")
end

end # module QMC

# --- To Run the Code ---
QMC.run_and_plot()