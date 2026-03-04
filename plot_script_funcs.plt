# --- Gnuplot Script ---

# 1. Set the output file. This will create a PNG image.
set terminal pngcairo enhanced font "Arial Black,14" size 640,480 linewidth 2
set output 'TEST.png'

# 2. Tell gnuplot the data files are comma-separated (CSV)
set datafile separator ","

# 3. Set plot titles and axis labels
set title "PA propagator, n=20, d=1, N=4" offset 0.1,-0.8
set xlabel "T"
set ylabel "C"

# 4. Set the axis ranges
set xrange [0:10]
set yrange [-5:1.5] 
#						          


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

# --- 6. Plotting Rules ---
# Use an even number of samples to avoid evaluating exactly at x=0
set samples 1000

# --- 7. The Plot Command ---
plot C_Tn(1/x, n, N)/n title "C_{Tn}(x)" with lines linewidth 2, \
     C_Hn(1/x, n, N)/n title "C_{Hn}(x)" with lines linewidth 2
