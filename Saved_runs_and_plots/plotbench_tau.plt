# 1. Define the Data Block
$Data << EOD
0.1	4.9066261349999
0.6	0.669453384000008
1.1	0.472590357999934
1.6	0.402417941000067
2.1	0.373585967000054
2.6	0.354249550000077
3.1	0.338402782999992
3.6	0.342090154000061
4.1	0.325258405999875
4.6	0.324200861000008
5.1	0.327131495000003
5.6	0.38405323500001
6.1	0.360968337000031
6.6	0.325030931000129
7.1	0.331784359000039
7.6	0.310739123000076
8.1	0.310747293999839
8.6	0.310131359000025
9.1	0.309533023999848
9.6	0.311683803999813
10.1	0.315097673000082
EOD

# 2. Define the Functions
g(x) = a*x**(-b) + c

# 3. Initialize Parameters
a = 0.17361387936111908
b = 1.4225030316326952
c = 0.31365774854470635

# 5. Plotting
set grid
set ytics 1
set xtics 2.5
set xrange [0:10]
set yrange [0:5]

set terminal postscript eps enhanced color font "Helvetica,28" linewidth 4 size 6,4.5
set output 'Bench_polyconv_tau.eps'

set encoding utf8
set xlabel "{{/Symbol t}"
set ylabel "t (s)"
set key top right Left

eq2 = sprintf("g(x) = %.2e * x^{-%.2f} + %.2f", a, b, c)

# Ensure no space follows the backslash below
plot $Data using 1:2 with points pt 7 ps 1.5 lc rgb "red" title "Computation", \
     g(x) with lines lw 2 lc rgb "blue" title eq2
