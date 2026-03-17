# 1. Define the Data Block (embedded data)

$Data << EOD
16	1250
32	3250
64	9500
96	17500
128	27500
160	38750
192	51000
224	64750
256	79500
288	95250
EOD



# 2. Define the Functions
# f(x) for Column 1 vs 2
f(x) = a1 * (x**1.5 - x)
g(x) = a2 * (x**1.5 + b2 * x)

# 3. Initialize Parameters (Crucial for power laws)
# Guessing quadratic (c=2) for the first plot
a1 = 21.36459

a2 = 20.4027889319
b2 = -1.06066017178

# 5. Plotting
set size ratio 1

set grid

set ytics 20000

set lmargin 1
set rmargin 1
set bmargin 2
set tmargin 2

set terminal postscript eps enhanced color font "Helvetica,20" linewidth 4
set output 'BitScaling_P.eps'




# --- Plot 1: n vs time ---
set title "Bit Precision Need Scaling"
set xlabel "n"
set ylabel "Simulated bits of precision needed"
set key top left Left reverse
# Create a label string with the found parameters
#eq1 = sprintf("f(x) =  %.2f * (x^{1.5} - x)" , a1)
eq1 = sprintf("g(x) =  %.2f * (x^{1.5} + %.2f x)" , a2, b2)
plot $Data using 1:2 with points pt 7 ps 1.5 title "Benchmark", \
	 g(x) with lines lw 2 title eq1




