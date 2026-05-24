# --- Gnuplot Script ---

# 1. Set output (increased the height slightly for the stacked layout)
set terminal postscript eps enhanced color font "Helvetica,24" linewidth 4 size 6,6
set output 'Mu_300_stacked.eps'

set datafile separator ","
set encoding utf8

# Begin multiplot mode
set multiplot

# ----------------------------------------------------
# --- MARGIN VARIABLES FOR PERFECT ALIGNMENT ---
# ----------------------------------------------------
# By defining these as variables, the left and right 
# edges of both graphs are guaranteed to align perfectly.
LM = 0.15   # Left margin
RM = 0.95   # Right margin
TM = 0.95   # Top margin of the canvas
MB = 0.35   # Middle boundary (where the two plots touch)
BM = 0.12   # Bottom margin of the canvas

# ==========================================
# --- TOP PLOT (Main Data) ---
# ==========================================
set lmargin at screen LM
set rmargin at screen RM
set tmargin at screen TM
set bmargin at screen MB

set title "PA propagator, n=300, d=2, T=0.01" offset 0,-0.5
set ylabel "Chemical potential" offset 1.5,0

# Hide X-axis labels and numbers for the top plot so they don't overlap the bottom plot
unset xlabel
set format x ""

set grid
unset key

# Labels
set label 1 "Exact" at 205, 19 font "Helvetica,20" textcolor rgb '#9400D3'
set label 2 "Thomas-Fermi" at 155, 16 font "Helvetica,20" textcolor rgb '#009E73'

# Line styles
set style line 3 lc rgb '#9400D3' lw 2
set style line 4 lc rgb '#009E73' lw 1

# Plot main data
plot 'chemical_potential_plot_data.csv' using 5:6 with lines ls 4, \
     ''                                 using 1:4 with lines ls 3


# ==========================================
# --- BOTTOM PLOT (Percentage Error) ---
# ==========================================
set lmargin at screen LM
set rmargin at screen RM
set tmargin at screen MB
set bmargin at screen BM

# Clear top labels and titles to prevent them from carrying over
unset label 1
unset label 2
unset title

# Restore X-axis labels for the bottom plot
set xlabel "n" offset 0,0.5
set format x "%h"  # Restores default number formatting on the x-axis

# Set Y-axis label for the smaller plot
set ylabel "% Error" offset 1.5,0

# Format: set ytics <start>, <step>, <end>
set ytics -5, 5, 15

# Plot the percentage error data
# tf_percentage_error_data.csv columns: n_TF, cp_pct_error, cv_pct_error
plot 'tf_percentage_error_data.csv' using 1:2 with lines lc rgb '#009E73' lw 0.8

# End multiplot
unset multiplot

# --- End of Script ---
