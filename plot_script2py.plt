# 1. Define the Data Block (embedded data)

$Data << EOD
2	0.746285851471298
4	0.535450504755717
8	0.276131080967832
16	0.097269877023798
32	0.027234565922072
64	0.00702256562581954
96	0.00313947909349634
128	0.00176959772787881
160	0.00113362351727373
192	0.000787646156575545
EOD



# 2. Define the Functions
# f(x) for Column 1 vs 2
f(x) = a1 / ( b1 + x**2)

# 3. Initialize Parameters (Crucial for power laws)
# Guessing quadratic (c=2) for the first plot
a1 = 28.81341
b1 = 36.22192

# 5. Plotting
set grid
set yrange [0:*]
set term qt size 1500,500 font "Arial Black,11.5" linewidth 2# Use 'qt' or 'wxt' or 'png'


# --- Plot 1: n vs time ---
set title "Relative Energy Error vs Number of Beads at Tau 15"
set xlabel "N"
set ylabel "Relative Energy Error"
set key top right Left reverse
# Create a label string with the found parameters
eq1 = sprintf("f(x) = %.2f /( %.2f + x^{2})", a1, b1)
plot $Data using 1:2 with points pt 7 ps 1.5 title "Benchmark", \
     f(x) with lines lw 2 title eq1


