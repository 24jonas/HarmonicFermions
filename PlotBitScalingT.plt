# 1. Define the Data Block (embedded data)

$Data << EOD
16	0.140266	0.0131
32	0.222928	0.0582
64	0.334747	0.3025
96	0.741987	0.8291
128	1.725691	1.7036
160	4.057761	3.0145
192	7.866986	4.8076
224	14.008861	6.9974
256	29.763399	9.8281
288	46.969451	13.2539
320	64.703063	17.6257
EOD



# 2. Define the Functions
# f(x) for Column 1 vs 2
f(x) = a1 + b1 * x**4
# g(x) for Column 1 vs 3
g(x) = a2 * x**b2

# 3. Initialize Parameters (Crucial for power laws)
# Guessing quadratic (c=2) for the first plot
a1 = 0.182267
b1 = 6.17707e-9

# Guessing quartic (c=4.5) for the second plot based on data shape
a2 = 0.0000137114
b2 = 2.42696

# 5. Plotting
set grid
set ytics 20
set yrange [0:*]
set terminal postscript eps enhanced color font "Helvetica,28" linewidth 4 size 6,4.5
set output 'BitScaling_T.eps'


# --- Plot 1: n vs time ---
set title "Time vs n"
set xlabel "n"
set ylabel "Time (s)"
set key top left Left reverse
# Create a label string with the found parameters
eq1 = sprintf("f(n) = %.3f + %.2e * n^4", a1, b1)
plot $Data using 1:2 with points pt 7 ps 2 title "Benchmark", \
	 f(x) with lines lw 2 title eq1

