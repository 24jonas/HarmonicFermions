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

$Data2 << EOD
10000	4822.92250474353
11000	4893.11186476536
12000	4220.86251341557
13000	4291.0518734374
14000	3990.02187777342
15000	3689.2929121051
16000	2831.273667322
17000	2715.71253509566
18000	2785.90189511749
19000	2484.87189945351
20000	1998.38602484765
21000	1883.28902938027
22000	1025.39472333378
23000	1095.2960378557
24000	608.803595443458
25000	678.992955465287
26000	377.661929805642
27000	0
28000	-224.328760938147
EOD



# 2. Define the Functions
# f(x) for Column 1 vs 2
f(x) = a1 * (x**1.5 - x)
# g(x) for Column 1 vs 3
g(x) = a2 + b2 * x

# 3. Initialize Parameters (Crucial for power laws)
# Guessing quadratic (c=2) for the first plot
a1 = 21.36459

# Guessing quartic (c=4.5) for the second plot based on data shape
a2 = 8060
b2 = -0.30103

# 5. Plotting
set grid

set ytics 1000
set xtics 4000
set xrange [10000:28000]
set yrange [-1000:5000]

#set lmargin 1
#set rmargin 1
#set bmargin 2
#set tmargin 2


set terminal postscript eps enhanced color font "Helvetica,28" linewidth 4 size 6,4.5
set output 'BitScaling_E.eps'



# --- Plot 2: n vs memory ---

set encoding utf8


set title "Δ_Z vs B for n=128"
set xlabel "B"
set ylabel "Δ_Z"
set key bottom left Left reverse
eq2 = sprintf("g(n) = 8060 - 0.301 * n")
plot $Data2 using 1:2 with points pt 7 ps 1.5 lc rgb "red" title "Computation", \
	 g(x) with lines lw 2 lc rgb "blue" title eq2


