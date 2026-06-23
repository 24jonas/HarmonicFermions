# 1. Define the Data Block (embedded data)

$Data << EOD
16	1250 1000 1250
32	3250 3000 3250
64	9500 9250 9500
96	17500 17250 17500
128	27500 27250 27500
160	38750 38500 38750
192	51000 50750 51000
224	64750 64500 64750
256	79500 79250 79500
288	95250 95000 95250
EOD



# 2. Define the Functions
# f(x) for Column 1 vs 2
g(x) = a2 * (x**1.5 + b2 * x)*15

# 3. Initialize Parameters (Crucial for power laws)
a2 = 1.3601859288
b2 = -1.06066017178

# 5. Plotting

set grid

set ytics 20000

#set lmargin 1
#set rmargin 1
#set bmargin 2
#set tmargin 2

set terminal postscript eps enhanced color font "Helvetica,26" linewidth 4 size 6,4.5
set output 'BitScaling_P.eps'




# --- Plot 1: n vs time ---
#set title "B vs n"
set xlabel "n"
set ylabel "B_2"
set key top left Left reverse
# Create a label string with the found parameters
#eq1 = sprintf("f(x) =  %.2f * (x^{1.5} - x)" , a1)
eq1 = sprintf("B_2(n,15)")
plot $Data using 1:2:3:4 with yerrorbars pt 7 ps 1.5 lw 5 title "Computation", \
	 g(x) with lines lw 2 title eq1




