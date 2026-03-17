# --- Gnuplot Script ---


$DataH << EOD
1.00000, 967.45816, .12521
1.50000, 956.38108, .07870
2.00000, 956.03086, .10159
2.50000, 960.82223, .12587
3.00000, 970.72304, .17025
3.50000, 984.98883, .21648
EOD

$DataT << EOD
1.00000, 938.27481, .18625
1.50000, 893.56622, .13478
2.00000, 849.58184, .13141
2.50000, 803.97907, .14226
3.00000, 757.36442, .15833
3.50000, 711.88985, .16724
EOD



# 1. Set the output file. This will create a PNG image.
set terminal postscript eps enhanced color font "Helvetica,22" linewidth 4
set output 'TEST.eps'

# 2. Tell gnuplot the data files are comma-separated (CSV)
set datafile separator ","

# 3. Set plot titles and axis labels
set title "PA propagator, n=100, d=2, N=[2,4,8,16]" offset 0.1,-0.8
set xlabel "tau"
set ylabel "E"

# 4. Set the axis ranges
set xrange [0:4]
set yrange [700:1000] 
#						          


set grid

# 5. Remove the legend (key)
unset key

# Labels for the Blue (T) curves
set label "T16" at 6.0, 935 font "Helvetica,16" textcolor rgb '#0000FF'
set label "T8" at 4.0, 925 font "Helvetica,16" textcolor rgb '#0000FF'
set label "T4" at 2.5, 915 font "Helvetica,16" textcolor rgb '#0000FF'
set label "T2" at 0.5, 905 font "Helvetica,16" textcolor rgb '#0000FF'

# Labels for the Red (H) curves
set label "H16" at 11.0, 952 font "Helvetica,16" textcolor rgb '#FF0000'
set label "H8" at 6.0, 952 font "Helvetica,16" textcolor rgb '#FF0000'
set label "H4" at 2.8, 954 font "Helvetica,16" textcolor rgb '#FF0000'
set label "H2" at 1.8, 965 font "Helvetica,16" textcolor rgb '#FF0000'



# 6. Define line styles with specific colors
# --- Shades of Red (for ham_PA_N1D1.csv) ---
set style line 3 lc rgb '#FF0000' # Red
set style line 4 lc rgb '#DC143C' # Crimson
set style line 5 lc rgb '#B22222' # Firebrick
set style line 6 lc rgb '#8B0000' # DarkRed

# --- Shades of Blue (for thermo_PA_N1D1.csv) ---
set style line 8 lc rgb '#1E90FF' # DodgerBlue
set style line 9 lc rgb '#4169E1' # RoyalBlue
set style line 10 lc rgb '#0000FF' # Blue
set style line 11 lc rgb '#0000CD' # MediumBlue
set style line 12 lc rgb '#00008B' # DarkBlue

# --- Dashed Black Line (for third file) ---
# dt 2 specifies a dash pattern.
set style line 13 lc rgb '#000000' dt 4

# 7. Plot the data (without title attributes)
plot 'data_comp_ham _N100 _D2 _PA.csv' using 1:2 with lines ls 6, \
     ''                 using 1:3 with lines ls 5, \
     ''                 using 1:4 with lines ls 4, \
	 ''					using 1:5 with lines ls 3, \
     'data_comp_thermo _N100 _D2 _PA.csv' using 1:2 with lines ls 12, \
     ''                   using 1:3 with lines ls 11, \
     ''                   using 1:4 with lines ls 10, \
	 ''					  using 1:5 with lines ls 9, \
     945 with lines ls 13, \
	 '$DataT' using 1:2:3 with yerrorbars pt 7 ps 1.5 lc rgb '#228B22', \
     '$DataH' using 1:2:3 with yerrorbars pt 7 ps 1.5 lc rgb '#00FF00'

# --- End of Script ---
