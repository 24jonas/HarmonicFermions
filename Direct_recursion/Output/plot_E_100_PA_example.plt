# --- Gnuplot Script ---

# 1. Set the output file. This will create a PNG image.
set terminal pngcairo enhanced color font "Helvetica,16" linewidth 2 size 600,450
set output 'E_100_PA.png'
#set output 'TEST.eps'

# 2. Tell gnuplot the data files are comma-separated (CSV)
set datafile separator ","

# 3. Set plot titles and axis labels
set encoding utf8
set title "PA propagator, n=100, d=2, N=[2,4,8,16]" offset 0.1,-0.8
set xlabel "{{/Symbol t}"
set ylabel "E"

# 4. Set the axis ranges
set xrange [0:10]
set yrange [900:990] 
#						          


set grid

set ytics 20

# 5. Remove the legend (key)
unset key

# Labels for the Blue (T) curves
set label "T16" at 6.0, 935 font "Helvetica,16" textcolor rgb '#0000FF'
set label "T8" at 4.0, 925 font "Helvetica,16" textcolor rgb '#0000FF'
set label "T4" at 2.5, 915 font "Helvetica,16" textcolor rgb '#0000FF'
set label "T2" at 0.5, 905 font "Helvetica,16" textcolor rgb '#0000FF'

# Labels for the Red (H) curves
set label "H16" at 11.0, 952 font "Helvetica,16" textcolor rgb '#FF0000'
set label "H8" at 6.3, 952 font "Helvetica,16" textcolor rgb '#FF0000'
set label "H4" at 3.1, 954 font "Helvetica,16" textcolor rgb '#FF0000'
set label "H2" at 2.2, 972 font "Helvetica,16" textcolor rgb '#FF0000'



# 6. Define line styles with specific colors
# --- Shades of Red (for ham_PA_N1D1.csv) ---
set style line 3 lc rgb '#FF0000' lw 2# Red
set style line 4 lc rgb '#FF0000' lw 2# Crimson
set style line 5 lc rgb '#FF0000' lw 2# Firebrick
set style line 6 lc rgb '#FF0000' lw 2# DarkRed

# --- Shades of Blue (for thermo_PA_N1D1.csv) ---
set style line 8 lc rgb '#0000FF' dt 5 lw 2# DodgerBlue
set style line 9 lc rgb '#0000FF' dt 5 lw 2 # RoyalBlue
set style line 10 lc rgb '#0000FF' dt 5 lw 2# Blue
set style line 11 lc rgb '#0000FF' dt 5 lw 2# MediumBlue
set style line 12 lc rgb '#0000FF' dt 5 lw 2# DarkBlue

# --- Dashed Black Line (for third file) ---
# dt 2 specifies a dash pattern.
set style line 13 lc rgb '#000000' dt 2 lw 2

# 7. Plot the data (without title attributes)
plot 'data_comp_ham _N100 _D2 _PA.csv' using 1:2 with lines ls 6, \
     ''                 using 1:3 with lines ls 5, \
     ''                 using 1:4 with lines ls 4, \
	 ''					using 1:5 with lines ls 3, \
     'data_comp_thermo _N100 _D2 _PA.csv' using 1:2 with lines ls 12, \
     ''                   using 1:3 with lines ls 11, \
     ''                   using 1:4 with lines ls 10, \
	 ''					  using 1:5 with lines ls 9, \
     945 with lines ls 13  # <--- This creates the horizontal dashed line

# --- End of Script ---
