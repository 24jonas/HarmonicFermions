# --- Gnuplot Script ---

# 1. Set the output file. This will create a PNG image.
set terminal pngcairo enhanced font "Arial Black,14" size 640,480 linewidth 2
set output 'all_curves_plot_3files.png'

# 2. Tell gnuplot the data files are comma-separated (CSV)
set datafile separator ","

# 3. Set plot titles and axis labels
set title "TI propagator, n=3, d=3" offset 0.1,-0.8
set xlabel "tau"
set ylabel "E"

# 4. Set the axis ranges
set xrange [0:20]
set yrange [6.45:6.6] # 0.45:0.6  1.95:2.1   4.45:4.6   12.45:12.6
#						          4.95:5.1   6.45:6.6


set grid

# 5. Remove the legend (key)
unset key

# 6. Define line styles with specific colors
# --- Shades of Red (for ham_PA_N1D1.csv) ---
set style line 1 lc rgb '#F08080' # LightCoral
set style line 2 lc rgb '#FF0000' # Red
set style line 3 lc rgb '#DC143C' # Crimson
set style line 4 lc rgb '#CD5C5C' # IndianRed
set style line 5 lc rgb '#B22222' # Firebrick
set style line 6 lc rgb '#8B0000' # DarkRed

# --- Shades of Blue (for thermo_PA_N1D1.csv) ---
set style line 7 lc rgb '#87CEFA' # LightSkyBlue
set style line 8 lc rgb '#1E90FF' # DodgerBlue
set style line 9 lc rgb '#4169E1' # RoyalBlue
set style line 10 lc rgb '#0000FF' # Blue
set style line 11 lc rgb '#0000CD' # MediumBlue
set style line 12 lc rgb '#00008B' # DarkBlue

# --- Dashed Black Line (for third file) ---
# dt 2 specifies a dash pattern.
set style line 13 lc rgb '#000000' dt 2

# 7. Plot the data (without title attributes)
plot 'data_ham_TI_N3D3.csv' using 1:2 with lines ls 1, \
     ''                 using 1:3 with lines ls 2, \
     ''                 using 1:4 with lines ls 3, \
     ''                 using 1:5 with lines ls 4, \
     ''                 using 1:6 with lines ls 5, \
     ''                 using 1:7 with lines ls 6, \
     'data_thermo_TI_N3D3.csv' using 1:2 with lines ls 7, \
     ''                   using 1:3 with lines ls 8, \
     ''                   using 1:4 with lines ls 9, \
     ''                   using 1:5 with lines ls 10, \
     ''                   using 1:6 with lines ls 11, \
     ''                   using 1:7 with lines ls 12, \
     'data_exact_N3D3.csv'   using 1:2 with lines ls 13 # <-- EDIT THIS FILENAME

# --- End of Script ---