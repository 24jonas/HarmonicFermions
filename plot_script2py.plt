# --- Gnuplot Script ---

# 1. Set the output file. This will create a PNG image.
set terminal pngcairo enhanced font "Arial Black,14" size 640,480 linewidth 2
set output 'all_curves_plot_3files.png'

# 2. Tell gnuplot the data files are comma-separated (CSV)
set datafile separator ","

# 3. Set plot titles and axis labels
set title "PA, TI, 4A propagator comparison in d=2" offset 0.1,-0.8
set xlabel "tau"
set ylabel "E - Exact Value"

# 4. Set the axis ranges
set xrange [0:20]
set yrange [-0.05:0.1] 


set grid

# 5. Remove the legend (key)
unset key

# 6. Define line styles with specific colors
set style line 1 lc rgb '#F08080' # LightCoral
set style line 2 lc rgb '#FF0000' # Red
set style line 3 lc rgb '#DC143C' # Crimson

set style line 4 lc rgb '#1E90FF' # DodgerBlue
set style line 5 lc rgb '#4169E1' # RoyalBlue
set style line 6 lc rgb '#0000FF' # Blue

set style line 7 lc rgb '#32CD32' # LimeGreen
set style line 8 lc rgb '#008000' # Green
set style line 9 lc rgb '#228B22' # ForestGreen

# --- Dashed Black Line (for third file) ---
# dt 2 specifies a dash pattern.
# set style line 13 lc rgb '#000000' dt 2

# 2 for d1, 3 for d2
# 12.5 for d1, 11 for d2 
# 50 for d1, 30 for d2

# 7. Plot the data (without title attributes)
plot 'data_comp_ham_N2_D2_PA.csv' using 1:($2 - 3) with lines ls 1, \
     'data_comp_ham_N5_D2_PA.csv' using 1:($2 - 11) with lines ls 2, \
	 'data_comp_ham_N10_D2_PA.csv' using 1:($2 - 30) with lines ls 3, \
     'data_comp_ham _N2 _D2 _TI.csv' using 1:($2 - 3) with lines ls 4, \
	 'data_comp_ham _N5 _D2 _TI.csv' using 1:($2 - 11) with lines ls 5, \
	 'data_comp_ham _N10 _D2 _TI.csv' using 1:($2 - 30) with lines ls 6, \
     'data_comp_ham_N2_D2_FA.csv' using 1:($2 - 3) with lines ls 7, \
	 'data_comp_ham_N5_D2_FA.csv' using 1:($2 - 11) with lines ls 8, \
	 'data_comp_ham_N10_D2_FA.csv' using 1:($2 - 30) with lines ls 9, \

# --- End of Script ---