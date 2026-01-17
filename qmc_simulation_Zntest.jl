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

using Statistics
using Printf

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


function get_ZnTau(bigfloat_precision::Int, propagator_choice::Symbol, tau_values::AbstractVector, bead_counts::AbstractVector, num_fermions::Int, dimensions::Int)
    setprecision(BigFloat, bigfloat_precision)

    # --- Propagator Definition ---
    p_funcs = propagators[propagator_choice]

    # --- Main Loop ---
    for N in bead_counts

        # 1. Pre-calculate b values (Thread-safe)
        b_values = Vector{BigFloat}(undef, length(tau_values))

        @threads for i in eachindex(tau_values)
            tau = tau_values[i]
            epsilon = BigFloat(tau) / BigFloat(N)
            
            # --- Standard Propagator (w=1) ---
            zeta_1 = p_funcs.zeta_1(epsilon)
            u = (zeta_1 >= 1) ? acosh(zeta_1) : BigFloat(0)
            b_values[i] = exp(-N * u)
        end

        num_points = length(b_values)
        Z_values = Vector{BigFloat}(undef, num_points)
        
        # --- Step A: Calculate Z for every point ---
        # We use a fresh cache for each point because 'b' is different for every tau.
        # Parallelizing this loop significantly speeds up high-precision calculation.
        Base.Threads.@threads for i in 1:num_points
            local_cache = Dict{Int, BigFloat}()
            Z_values[i] = Z_recursive!(num_fermions, dimensions, b_values[i], local_cache)
        end

        return Z_values
    end

    
end


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

    test_precision = 34000 #starting point
    precision_step = 1000
    err_threshold = 0.00000001

    ref_Zn = get_ZnTau(bigfloat_precision, propagator_choice, tau_values, bead_counts, num_fermions, dimensions)
    test_Zn = get_ZnTau(test_precision, propagator_choice, tau_values, bead_counts, num_fermions, dimensions)
    err = abs.((ref_Zn .- test_Zn) ./ ref_Zn)
    max_err = maximum(err)

    while isnan(max_err) || max_err > err_threshold
        test_precision = test_precision + precision_step
        test_Zn = get_ZnTau(test_precision, propagator_choice, tau_values, bead_counts, num_fermions, dimensions)
        err = abs.((ref_Zn .- test_Zn) ./ ref_Zn)
        max_err = maximum(err)

        print(" ---- ")
        @printf "%.4e\n" max_err
        print(" <--> ")
        print(test_precision)
        print(" ---- ")
    end

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