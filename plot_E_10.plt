# --- Gnuplot Script ---

# 1. Set the output file. This will create a PNG image.
set terminal postscript eps enhanced color font "Helvetica,22" linewidth 4
set output 'E_forC.eps'

# 2. Tell gnuplot the data files are comma-separated (CSV)
set datafile separator ","

# 3. Set plot titles and axis labels
set title "PA propagator, n=10, d=1, N=[2,20,200]" offset 0.1,-0.8
set xlabel "tau"
set ylabel "E"

# 4. Set the axis ranges
set xrange [0:4]
set yrange [45:55] 
#						          


set grid

# 5. Remove the legend (key)
unset key

# Labels for the Blue (T) curves
set label "T200" at 3.0, 50.3 font "Helvetica,16" textcolor rgb '#0000FF'
set label "T20" at 3, 49.5 font "Helvetica,16" textcolor rgb '#0000FF'
set label "T2" at 1, 48 font "Helvetica,16" textcolor rgb '#0000FF'

# Labels for the Red (H) curves
set label "H200" at 3.7, 50.3 font "Helvetica,16" textcolor rgb '#FF0000'
set label "H20" at 3.35, 50.3 font "Helvetica,16" textcolor rgb '#FF0000'
set label "H2" at 2, 51 font "Helvetica,16" textcolor rgb '#FF0000'



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
plot 'data_comp_ham _N10 _D1 _PA.csv' using 1:2 with lines ls 6, \
     ''                 using 1:3 with lines ls 5, \
     ''                 using 1:4 with lines ls 4, \
	 ''					using 1:5 with lines ls 3, \
     'data_comp_thermo _N10 _D1 _PA.csv' using 1:2 with lines ls 12, \
     ''                   using 1:3 with lines ls 11, \
     ''                   using 1:4 with lines ls 10, \
	 ''					  using 1:5 with lines ls 9, \
     50 with lines ls 13  # <--- This creates the horizontal dashed line

# --- End of Script ---
