using Plots
using CSV
using DataFrames
using JSON

# The function now takes the unique run_id (the timestamp string) as input
function create_plot(run_id::String)
    println("Generating plot for run: '$run_id'")

    # --- Construct filenames using the run_id ---
    param_file = "parameters_$(run_id).json"
    data_file = "data_$(run_id).csv"

    if !isfile(param_file) || !isfile(data_file)
        println("Error: Could not find required files for run_id '$run_id'.")
        return
    end

    # --- Load data and parameters (logic is the same) ---
    params = JSON.parsefile(param_file)
    df = CSV.read(data_file, DataFrame)

    num_fermions = params["num_fermions"]
    dimensions = params["dimensions"]
    
    # --- Plotting Logic (unchanged) ---
    plt = plot(
        title = "Energy vs. Tau for n=$num_fermions, d=$dimensions",
        xlabel = "τ (Imaginary Time)",
        ylabel = "Energy E",
        legend = :topright,
        grid = true
    )

    tau_values = df.tau
    for col_name in names(df)[2:end]
        label_text = replace(col_name, "_" => " = ") * " beads"
        plot!(plt, tau_values, df[!, col_name], label=label_text)
    end

    # --- Display and Save the Plot ---
    display(plt)
    output_plot_filename = "plot_$(run_id).png"
    savefig(plt, output_plot_filename)
    println("Plot saved to '$output_plot_filename'.")
end

function compare_runs(run_id1::String, run_id2::String)
    println("Generating comparison plot for runs:")
    println("   Run 1: $run_id1")
    println("   Run 2: $run_id2")

    # --- Load Data and Parameters for Run 1 ---
    param_file1 = "parameters_$(run_id1).json"
    data_file1 = "data_$(run_id1).csv"
    if !isfile(param_file1) || !isfile(data_file1)
        println("Error: Could not find required files for run_id '$run_id1'.")
        return
    end
    params1 = JSON.parsefile(param_file1)
    df1 = CSV.read(data_file1, DataFrame)

    # --- Load Data and Parameters for Run 2 ---
    param_file2 = "parameters_$(run_id2).json"
    data_file2 = "data_$(run_id2).csv"
    if !isfile(param_file2) || !isfile(data_file2)
        println("Error: Could not find required files for run_id '$run_id2'.")
        return
    end
    params2 = JSON.parsefile(param_file2)
    df2 = CSV.read(data_file2, DataFrame)

    # --- 1. Define Color Palettes for Each Run ---
    # We'll generate a gradient of blues and reds. The number of colors
    # should match the number of lines we need to plot for each run.
    # We add 2 and slice [3:end] to avoid the lightest colors, which can be hard to see.
    num_lines1 = length(names(df1)) - 1
    num_lines2 = length(names(df2)) - 1
    blues = palette(:blues, num_lines1 + 2)[3:end]
    reds = palette(:reds, num_lines2 + 2)[3:end]

    # --- Initialize the Plot ---
    plt = plot(
        title = "Energy vs Tau for n=3, d=2",
        xlabel = "τ",
        ylabel = "Energy E",
        legend = :topright,
        grid = true,
        size = (1200, 800),  # <-- Makes the image 1200px wide and 800px tall
        dpi = 300,
        ylim = (2, 14)
    )

    # --- Plot Data for Run 1 (shades of blue, solid lines) ---
    tau_values1 = df1.tau
    # --- 2. Loop with an index using enumerate ---
    for (i, col_name) in enumerate(names(df1)[2:end])
        label_text = "$(run_id1): " * replace(col_name, "_" => " = ")
        # --- 3. Assign the specific color and linestyle ---
        plot!(plt, tau_values1, df1[!, col_name], label=label_text, linestyle=:solid, color=blues[i])
    end

    # --- Plot Data for Run 2 (shades of red, dashed lines) ---
    tau_values2 = df2.tau
    for (i, col_name) in enumerate(names(df2)[2:end])
        label_text = "$(run_id2): " * replace(col_name, "_" => " = ")
        plot!(plt, tau_values2, df2[!, col_name], label=label_text, linestyle=:dash, color=reds[i])
    end

    # --- Display and Save the Plot ---
    display(plt)
    
    output_plot_filename = "comparison_$(run_id1)_vs_$(run_id2).png"
    savefig(plt, output_plot_filename)
    println("✅ Comparison plot saved to '$output_plot_filename'.")
end




# Filename: plot_results.jl

# ... (using statements: Plots, CSV, DataFrames, JSON) ...

# This function is called by `plot_latest_run`.
# It does the actual work of loading files and creating the plot for a given run_id.
function create_plot(run_id::String)
    println("🎨 Generating plot for run: '$run_id'")

    # --- Construct filenames using the run_id ---
    param_file = "parameters_$(run_id).json"
    data_file = "data_$(run_id).csv"

    if !isfile(param_file) || !isfile(data_file)
        println("Error: Could not find required files for run_id '$run_id'.")
        return
    end

    # --- Load data and parameters ---
    params = JSON.parsefile(param_file)
    df = CSV.read(data_file, DataFrame)
    
    # --- Plotting Logic ---
    plt = plot(
        title = "Energy vs. Tau (Run: $(run_id))",
        xlabel = "τ (Imaginary Time)",
        ylabel = "Energy E",
        legend = :topright,
        grid = true,
        size = (1200, 800), # Example size
        dpi = 300           # Example resolution
    )
    
    tau_values = df.tau
    for col_name in names(df)[2:end]
        label_text = replace(col_name, "_" => " = ") * " beads"
        plot!(plt, tau_values, df[!, col_name], label=label_text)
    end
    
    # --- Display and Save the Plot ---
    display(plt)
    output_plot_filename = "plot_$(run_id).png"
    savefig(plt, output_plot_filename)
    println("✅ Plot saved to '$output_plot_filename'.")
end