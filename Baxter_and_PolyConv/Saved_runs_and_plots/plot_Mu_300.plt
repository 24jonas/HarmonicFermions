# --- Gnuplot Script ---

# 1. Set the output file. This will create a PNG image.
set terminal postscript eps enhanced color font "Helvetica,28" linewidth 4 size 6,4.5
set output 'Mu_300.eps'

# 2. Tell gnuplot the data files are comma-separated (CSV)
set datafile separator ","

# 3. Set plot titles and axis labels
set encoding utf8
set title "PA propagator, n=300, d=2, T=0.01" offset 0.1,-0.8
set xlabel "n"
set ylabel "Chemical potential"



# 4. Set the axis ranges
#set xrange [0:15]
#set yrange [0.95:1.05] 
						          


set grid
#set ytics Eexact*0.02
#set format y "%.0f"

# 5. Remove the legend (key)
unset key

# Labels for the Blue (T) curves
set label "Exact" at 205, 19 font "Helvetica,20" textcolor rgb '#000000'
set label "Thomas-Fermi" at 155, 16 font "Helvetica,20" textcolor rgb '#FF0000'



# 6. Define line styles with specific colors
# --- Shades of Red (for ham_PA_N1D1.csv) ---
set style line 3 lc rgb '#000000' lw 2
set style line 4 lc rgb '#FF0000' lw 1

# --- Dashed Black Line (for third file) ---
# dt 2 specifies a dash pattern.
set style line 13 lc rgb '#000000' dt 2 lw 2

# 7. Plot the data (without title attributes)
plot 'chemical_potential_plot_data.csv' using 1:4 with lines ls 3, \
     ''                 using 5:6 with lines ls 4, \

# --- End of Script ---
