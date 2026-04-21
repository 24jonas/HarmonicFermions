# 1. Define the Data Block (embedded data)

$Data << EOD
0	2.92000004264992E-06
250	0.0257673700000396
500	0.0864605230000279
750	0.168187095999997
1000	0.282075028999998
1250	0.411116303999961
1500	0.559229519000041
1750	0.747509654999931
2000	0.921539640999981
2250	1.16716383899984
2500	1.40604885800008
2750	1.65859313300007
3000	1.9728496329999
3250	2.32110423299991
3500	2.61606040400011
3750	2.9559576040001
4000	3.38519906900001
4250	3.65499519500008
4500	4.11746982499994
4750	4.64074525699994
5000	5.66782104499998
EOD



# 2. Define the Functions
# g(x) for Column 1 vs 3
g(x) = a*x^b + c

# 3. Initialize Parameters (Crucial for power laws)
# Guessing quartic (c=4.5) for the second plot based on data shape
a = 1.8325333423707465e-07
b = 2.0156759585326016
c = 0.06574591211478839

# 5. Plotting
set grid

set ytics 1
set xtics 1000
set xrange [0:5000]
set yrange [0:6]

#set lmargin 1
#set rmargin 1
#set bmargin 2
#set tmargin 2


set terminal postscript eps enhanced color font "Helvetica,28" linewidth 4 size 6,4.5
set output 'Bench_polyconv.eps'



# --- Plot 2: n vs memory ---

set encoding utf8


#set title "Δ_Z vs B for n=128"
set xlabel "{/Symbol b}"
set ylabel "Δ_Z^2"
set key bottom left Left reverse
eq2 = sprintf("g(x) = %.2e * x^{%.2f} + %.2f", a, b, c)
plot $Data2 using 1:2 with points pt 7 ps 1.5 lc rgb "red" title "Computation", \
	 g(x) with lines lw 2 lc rgb "blue" title eq2
