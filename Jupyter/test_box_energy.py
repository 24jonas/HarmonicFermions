import math

def logaddexp(a, b):
    if a == -math.inf: return b
    if b == -math.inf: return a
    if a > b: return a + math.log1p(math.exp(b - a))
    else: return b + math.log1p(math.exp(a - b))

def degeneracy_box(d, k):
    if d == 1:
        root = int(math.sqrt(k))
        if root * root == k and root > 0: return 1
        return 0
    elif d == 2:
        count = 0
        max_n = int(math.sqrt(k))
        for nx in range(1, max_n + 1):
            ny_sq = k - nx**2
            if ny_sq > 0:
                ny = int(math.sqrt(ny_sq))
                if ny**2 == ny_sq: count += 1
        return count
    return 0

def log_binom(n, k):
    if k < 0 or k > n: return -math.inf
    return math.lgamma(n + 1) - math.lgamma(k + 1) - math.lgamma(n - k + 1)

def fermion_logZ_numeric(tau, n, d, L, k_start=1):
    logZ = [-math.inf] * (n + 1)
    logZ[0] = 0.0
    
    def convolve_shell(logZ, g, logw, n):
        mmax = min(g, n)
        logC = [0.0] * (mmax + 1)
        for m in range(1, mmax + 1):
            logC[m] = log_binom(g, m) + m * logw
        old = logZ[:]
        for r in range(n + 1):
            s = -math.inf
            upper = min(r, mmax)
            for m in range(upper + 1):
                if old[r - m] != -math.inf:
                    s = logaddexp(s, logC[m] + old[r - m])
            logZ[r] = s

    cumulative = 0
    k = k_start
    while cumulative < n or k < k_start + 1000:
        g = degeneracy_box(d, k)
        if g > 0:
            E_k = (math.pi**2 / (2 * L**2)) * k
            logw = -tau * E_k
            convolve_shell(logZ, g, logw, n)
            cumulative += g
        k += 1
    return logZ[n]

v1 = fermion_logZ_numeric(2.0, 100, 2, 2.0)
v2 = fermion_logZ_numeric(2.5, 100, 2, 2.0)
print("E_GS for L=2.0:", (v1 - v2) / 0.5)

v3 = fermion_logZ_numeric(2.0, 100, 2, 4.0)
v4 = fermion_logZ_numeric(2.5, 100, 2, 4.0)
print("E_GS for L=4.0:", (v3 - v4) / 0.5)
