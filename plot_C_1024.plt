# --- Gnuplot Script ---

# 1. Set the output file. This will create a PNG image.
set terminal postscript eps enhanced color font "Helvetica,22" linewidth 4
set output 'C_1024_PA.eps'

# 2. Tell gnuplot the data files are comma-separated (CSV)
set datafile separator ","

# 3. Set plot titles and axis labels
set title "PA propagator, n=1024, d=2, N=[2,4,8,16], w=1/2" offset 0.1,-0.8
set xlabel "T"
set ylabel "C/n"

# 4. Set the axis ranges
set xrange [0:2]
set yrange [-0.5:1.5] 
#						          


set grid

# 5. Remove the legend (key)
unset key

# Labels for the Blue (T) curves
set label "T2" at 1.0, 0.65 font "Helvetica,16" textcolor rgb '#0000FF'
set label "T4" at 0.65, 0.5 font "Helvetica,16" textcolor rgb '#0000FF'
set label "T8" at 0.405, 0.42 font "Helvetica,16" textcolor rgb '#0000FF'
set label "T16" at 0.275, 0.3 font "Helvetica,16" textcolor rgb '#0000FF'

# Labels for the Red (H) curves
set label "H16" at 0.04, 0.1 font "Helvetica,16" textcolor rgb '#FF0000'
set label "H8" at 0.18, 0.17 font "Helvetica,16" textcolor rgb '#FF0000'
set label "H4" at 0.23, -0.25 font "Helvetica,16" textcolor rgb '#FF0000'
set label "H2" at 0.35, -0.4 font "Helvetica,16" textcolor rgb '#FF0000'



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
plot 'data_Cv_comp_ham _N1024 _D2 _PA.csv' using 1:2 with lines ls 6, \
     ''                 using 1:3 with lines ls 5, \
     ''                 using 1:4 with lines ls 4, \
	 ''					using 1:5 with lines ls 3, \
     'data_Cv_comp_thermo _N1024 _D2 _PA.csv' using 1:2 with lines ls 12, \
     ''                   using 1:3 with lines ls 11, \
     ''                   using 1:4 with lines ls 10, \
	 ''					  using 1:5 with lines ls 9, \
     # 945 with lines ls 13  # <--- This creates the horizontal dashed line

# --- End of Script ---
