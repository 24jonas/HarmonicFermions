# 1. Define the Data Block (embedded data)

$Data << EOD
16	0.25	0.09
32	0.28	0.28
64	0.93	1.33
96	4.1		3.43
128	14.96	6.99
160	37.6	12.08
192	81.42	19.06
224	165.12	28.19
256	305.38	39.49
288	509.26	53.18
EOD



# 2. Define the Functions
# f(x) for Column 1 vs 2
f(x) = a1 + b1 * x**c1
# g(x) for Column 1 vs 3
g(x) = a2 + b2 * x**c2

# 3. Initialize Parameters (Crucial for power laws)
# Guessing quadratic (c=2) for the first plot
a1 = 0.249067
b1 = 7.2922e-9
c1 = 4.40713

# Guessing quartic (c=4.5) for the second plot based on data shape
a2 = 0.0496673
b2 = 0.0000424963
c2 = 2.47625

# 5. Plotting
set grid
set yrange [0:*]
set term wxt size 1500,500 font "Arial Black,11.5" linewidth 2# Use 'qt' or 'wxt' or 'png'

set multiplot layout 1,2

    # --- Plot 1: n vs time ---
    set title "Computation Time Scaling"
    set xlabel "n"
    set ylabel "Time [s]"
    set key top left Left reverse
    # Create a label string with the found parameters
    eq1 = sprintf("f(x) = %.2f + %.2e * x^{%.2f}", a1, b1, c1)
    plot $Data using 1:2 with points pt 7 ps 1.5 title "Benchmark", \
         f(x) with lines lw 2 title eq1

    # --- Plot 2: n vs memory ---
	set title "Memory Usage Scaling"
    set xlabel "n"
	set ylabel "Peak [MB] usage by n-size dependant arrays"
    set key top left Left reverse
    eq2 = sprintf("g(x) = %.2f + %.2e * x^{%.2f}", a2, b2, c2)
    plot $Data using 1:3 with points pt 7 ps 1.5 lc rgb "red" title "Benchmark", \
         g(x) with lines lw 2 lc rgb "blue" title eq2

unset multiplot