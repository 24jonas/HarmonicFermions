# --- Gnuplot Script ---

# 1. Set the output file. This will create a PNG image.
set terminal postscript eps enhanced color font "Helvetica,28" linewidth 4 size 6,4.5
set output 'C_analytical.eps'

# 2. Tell gnuplot the data files are comma-separated (CSV)
set datafile separator ","

# 3. Set plot titles and axis labels
set title "PA propagator, n=10, d=1, N=[2,20,200]" offset 0.1,-0.8
set xlabel "T"
set ylabel "C/n"

# 4. Set the axis ranges
set xrange [0:0.8]
set yrange [-0.5:0.8] 
#						          

# Labels for the Blue (T) curves
set label "T200" at 0.1, 0.15 font "Helvetica,20" textcolor rgb '#0000FF'
set label "T20" at 0.25, 0.35 font "Helvetica,20" textcolor rgb '#0000FF'
set label "T2" at 1, 0.7 font "Helvetica,20" textcolor rgb '#0000FF'

# Labels for the Red (H) curves
set label "H200" at 0.03, -0.06 font "Helvetica,20" textcolor rgb '#FF0000'
set label "H20" at 0.12, -0.3 font "Helvetica,20" textcolor rgb '#FF0000'
set label "H2" at 0.82, 0.1 font "Helvetica,20" textcolor rgb '#FF0000'

# Label for reference
set label "Exact" at 0.4, 0.0 font "Helvetica,20" textcolor rgb '#000000'

set grid

# 5. Remove the legend (key)
# set key top right box  # Puts the legend in a box at the top right
unset key

# 1. Define E_n(x)
E(x, n, N) = n * (n - 1) / 4.0 + \
             0.5 * sum [i=1:n] (i * (1.0 / tanh(i * N * asinh(x / (2.0 * N)))))

# 2. Define C_Tn(x)
C_Tn(x, n, N) = x**2 * ( \
    (x / (4.0 * N**2 * (1.0 + x**2 / (4.0 * N**2))**1.5)) * E(x, n, N) + \
    (1.0 / (4.0 * (1.0 + x**2 / (4.0 * N**2)))) * \
    sum [i=1:n] (i**2 * (1.0 / sinh(i * N * asinh(x / (2.0 * N))))**2) \
)

# 3. Define C_Hn(x)
C_Hn(x, n, N) = x**2 * ( \
    ((8.0 * N**2 + x**2) / (8.0 * (4.0 * N**2 + x**2))) * \
    sum [i=1:n] (i**2 * (1.0 / sinh(i * N * asinh(x / (2.0 * N))))**2) - \
    (x**3 / (32.0 * N**4 * (1.0 + x**2 / (4.0 * N**2))**1.5)) * E(x, n, N) \
)


# --- 4. Define Parameters ---
# You can easily change these values to see how the curves behave
N = 5.0
n = 10


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

set style line 13 lc rgb '#000000' # ?

# --- 6. Plotting Rules ---
# Use an even number of samples to avoid evaluating exactly at x=0
set samples 1000

# --- 7. The Plot Command ---
plot C_Tn(1/x, n, 2)/n title "C_{Tn}(x)" with lines ls 8, \
     C_Tn(1/x, n, 20)/n title "C_{Tn}(x)" with lines ls 9, \
     C_Tn(1/x, n, 200)/n title "C_{Tn}(x)" with lines ls 10, \
     C_Hn(1/x, n, 2)/n title "C_{Hn}(x)" with lines ls 3, \
     C_Hn(1/x, n, 20)/n title "C_{Hn}(x)" with lines ls 4, \
     C_Hn(1/x, n, 200)/n title "C_{Hn}(x)" with lines ls 5, \
     C_Hn(1/x, n, 2000)/n title "C_{Hn}(x)" with lines ls 13
