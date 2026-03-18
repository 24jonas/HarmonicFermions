# --- Gnuplot Script ---

# 1. Set the output file. This will create a PNG image.
set terminal postscript eps enhanced color font "Helvetica,28" linewidth 4 size 6,6
set output 'E_100_4A.eps'
#set output 'TEST.eps'

# 2. Tell gnuplot the data files are comma-separated (CSV)
set datafile separator ","

# 3. Set plot titles and axis labels
set encoding utf8
set title "4A propagator, n=100, d=2, N=[2,4,8,16]" offset 0.1,-0.8
set xlabel "{{/Symbol t} ((\U+210F {/Symbol w})^{-1})"
set ylabel "E ( \U+210F {/Symbol w})"

# 4. Set the axis ranges
set xrange [0:15]
set yrange [900:990] 
#						          

set size ratio 1

set grid

set ytics 20

# 5. Remove the legend (key)
unset key

# Labels for the Blue (T) curves
set label "T16" at 13.7, 942 font "Helvetica,20" textcolor rgb '#0000FF'
set label "T8" at 12.0, 936 font "Helvetica,20" textcolor rgb '#0000FF'
set label "T4" at 10.0, 922 font "Helvetica,20" textcolor rgb '#0000FF'
set label "T2" at 6.0, 910 font "Helvetica,20" textcolor rgb '#0000FF'

# Labels for the Red (H) curves
#set label "H16" at 11.0, 952 font "Helvetica,20" textcolor rgb '#FF0000'
set label "H8" at 12.7, 942 font "Helvetica,20" textcolor rgb '#FF0000'
set label "H4" at 10.0, 948 font "Helvetica,20" textcolor rgb '#FF0000'
set label "H2" at 6.9, 952 font "Helvetica,20" textcolor rgb '#FF0000'



# 6. Define line styles with specific colors
# --- Shades of Red (for ham_PA_N1D1.csv) ---
set style line 3 lc rgb '#FF0000' lw 2# Red
set style line 4 lc rgb '#DC143C' lw 2# Crimson
set style line 5 lc rgb '#B22222' lw 2# Firebrick
set style line 6 lc rgb '#8B0000' lw 2# DarkRed

# --- Shades of Blue (for thermo_PA_N1D1.csv) ---
set style line 8 lc rgb '#1E90FF' dt 5 lw 2# DodgerBlue
set style line 9 lc rgb '#4169E1' dt 5 lw 2# RoyalBlue
set style line 10 lc rgb '#0000FF' dt 5 lw 2# Blue
set style line 11 lc rgb '#0000CD' dt 5 lw 2# MediumBlue
set style line 12 lc rgb '#00008B' dt 5 lw 2# DarkBlue

# --- Dashed Black Line (for third file) ---
# dt 2 specifies a dash pattern.
set style line 13 lc rgb '#000000' dt 2 lw 2

# 7. Plot the data (without title attributes)
plot 'data_comp_ham _N100 _D2 _FA.csv' using 1:2 with lines ls 6, \
     ''                 using 1:3 with lines ls 5, \
     ''                 using 1:4 with lines ls 4, \
	 ''					using 1:5 with lines ls 3, \
     'data_comp_thermo _N100 _D2 _FA.csv' using 1:2 with lines ls 12, \
     ''                   using 1:3 with lines ls 11, \
     ''                   using 1:4 with lines ls 10, \
	 ''					  using 1:5 with lines ls 9, \
     945 with lines ls 13  # <--- This creates the horizontal dashed line

# --- End of Script ---
