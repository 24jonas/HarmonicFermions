# 1. Define the Data Block (embedded data)

$Data << EOD
2.66666666666667	1.03904993539877
1.6	0.840401519555659
1.14285714285714	0.691149420032646
0.888888888888889	0.564152933864044
0.727272727272727	0.460476280572617
0.615384615384615	0.378371230637064
0.533333333333333	0.313925310480883
0.470588235294118	0.263044052491992
EOD

term(l, x) = (x == 0) ? 1.0 : ((l * x / 2.0)**2) / (sinh(l * x / 2.0)**2)
C4(x) = term(1, x) + term(2, x) + term(3, x) + term(4, x)

# 5. Plotting
set grid
set xrange [0:2.75]
set yrange [0:*]
set term qt size 800,600 font "Arial Black,11.5" linewidth 2# Use 'qt' or 'wxt' or 'png'


# --- Plot 1: n vs time ---
set title "Heat Capacity, PA, n=4, d=1, N=4*tau"
set xlabel "T"
set ylabel "(C_V)/n"
set key top left Left reverse
# Create a label string with the found parameters
plot $Data using 1:2 with points pt 7 ps 1.5 title "Benchmark", \
     C4(1/x)/4 with lines lw 2 title "E"


