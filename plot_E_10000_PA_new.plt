# --- Gnuplot Script ---

# 1. Set the output file. This will create a PNG image.
set terminal postscript eps enhanced color font "Helvetica,28" linewidth 4 size 6,4.5
set output 'E_10000_PA.eps'

# 2. Tell gnuplot the data files are comma-separated (CSV)
set datafile separator ","

# 3. Set plot titles and axis labels
set encoding utf8
set title "n=10,000, d=2, N=[2,4,8,16], {/Symbol w}=0.5" offset 0.1,-0.8
set xlabel "{{/Symbol t}"
set ylabel "E / 10^5"

Eexact = 471410.5
scale = 0.0008

# 4. Set the axis ranges
set xrange [0:10]
set yrange [Eexact*0.998:Eexact*1.002]
						          


set grid
set ytics Eexact*scale
set format y "%.3t"

# 5. Remove the legend (key)
unset key

# Labels for the Blue (T) curves
set label "T16" at 2.4, Eexact-0.4*Eexact*scale font "Helvetica,20" textcolor rgb '#0000FF'
set label "T8" at 2.1, Eexact-1.6*Eexact*scale font "Helvetica,20" textcolor rgb '#0000FF'
set label "T4" at 1.2, Eexact-1.8*Eexact*scale font "Helvetica,20" textcolor rgb '#0000FF'
set label "T2" at 0.2, Eexact-2.3*Eexact*scale font "Helvetica,20" textcolor rgb '#0000FF'

# Labels for the Red (H) curves
set label "H16" at 9, Eexact+0.3*Eexact*scale font "Helvetica,20" textcolor rgb '#FF0000'
set label "H8" at 5, Eexact+0.4*Eexact*scale font "Helvetica,20" textcolor rgb '#FF0000'
set label "H4" at 2.1, Eexact+0.6*Eexact*scale font "Helvetica,20" textcolor rgb '#FF0000'
set label "H2" at 1.5, Eexact+1.5*Eexact*scale font "Helvetica,20" textcolor rgb '#FF0000'



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
     ''                 using 1:8 with lines ls 8, \
     ''                   using 1:3 with lines ls 3, \
     ''                   using 1:5 with lines ls 3, \
	 ''					  using 1:7 with lines ls 3, \
	 ''					  using 1:9 with lines ls 3, \
     Eexact with lines ls 13  # <--- This creates the horizontal dashed line

# --- End of Script ---
