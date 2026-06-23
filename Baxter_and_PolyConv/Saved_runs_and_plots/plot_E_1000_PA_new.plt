# --- Gnuplot Script ---

# 1. Set the output file. This will create a PNG image.
set terminal postscript eps enhanced color font "Helvetica,28" linewidth 4 size 6,4.5
set output 'E_10000_PA.eps'

# 2. Tell gnuplot the data files are comma-separated (CSV)
set datafile separator ","

# 3. Set plot titles and axis labels
set encoding utf8
set title "PA propagator, n=10,000, d=2, N=[2,4,8]" offset 0.1,-0.8
set xlabel "{{/Symbol t}"
set ylabel "E"

Eexact = 942820.0
scale = 0.004

# 4. Set the axis ranges
set xrange [0:2.5]
set yrange [Eexact*0.99:Eexact*1.01]
						          


set grid
set ytics Eexact*0.004
set format y "%.3t · 10^{%T}"

# 5. Remove the legend (key)
unset key

# Labels for the Blue (T) curves
#set label "T16" at 6.0, Eexact-0.7*Eexact*0.02 font "Helvetica,20" textcolor rgb '#0000FF'
set label "T8" at 1.8, Eexact-1.3*Eexact*scale font "Helvetica,20" textcolor rgb '#0000FF'
set label "T4" at 1.1, Eexact-1.8*Eexact*scale font "Helvetica,20" textcolor rgb '#0000FF'
set label "T2" at 0.6, Eexact-2.3*Eexact*scale font "Helvetica,20" textcolor rgb '#0000FF'

# Labels for the Red (H) curves
#set label "H16" at 11.0, Eexact+0.4*Eexact*0.02 font "Helvetica,20" textcolor rgb '#FF0000'
set label "H8" at 2.3, Eexact-0.2*Eexact*scale font "Helvetica,20" textcolor rgb '#FF0000'
set label "H4" at 2.1, Eexact+0.4*Eexact*scale font "Helvetica,20" textcolor rgb '#FF0000'
set label "H2" at 1.7, Eexact+1.5*Eexact*scale font "Helvetica,20" textcolor rgb '#FF0000'



# 6. Define line styles with specific colors
# --- Shades of Red (for ham_PA_N1D1.csv) ---
set style line 3 lc rgb '#FF0000' lw 1.5# Red


# --- Shades of Blue (for thermo_PA_N1D1.csv) ---
set style line 8 lc rgb '#0000FF' dt 5 lw 1.5# DodgerBlue


# --- Dashed Black Line (for third file) ---
# dt 2 specifies a dash pattern.
set style line 13 lc rgb '#000000' dt 2 lw 2

# 7. Plot the data (without title attributes)
plot 'extracted_n10000_all_beads.csv' using 1:2 with lines ls 8, \
     ''                 using 1:4 with lines ls 8, \
     ''                 using 1:6 with lines ls 8, \
     ''                   using 1:3 with lines ls 3, \
     ''                   using 1:5 with lines ls 3, \
	 ''					  using 1:7 with lines ls 3, \
     Eexact with lines ls 13  # <--- This creates the horizontal dashed line

# --- End of Script ---
