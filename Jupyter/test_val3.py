import math
# using the exact logic but with level_degeneracy
def logaddexp(a, b):
    if a == -math.inf: return b
    if b == -math.inf: return a
    if a > b: return a + math.log1p(math.exp(b - a))
    else: return b + math.log1p(math.exp(a - b))
def log_binom(n, k):
    if k < 0 or k > n: return -math.inf
    return math.lgamma(n + 1) - math.lgamma(k + 1) - math.lgamma(n - k + 1)

def fermion_logZ_numeric(tau, N, n, d=None, max_shell=None, tol=1e-4, consecutive_small=8, safety_cap=100000, return_all=False, gk_fn=None, Ek_fn=None, logwk_fn=None, k_start=0):
    epsilon = tau / N
    zeta = 1.0 + 0.5 * epsilon * epsilon
    u = math.log(zeta + math.sqrt(zeta * zeta - 1.0))
    b = math.exp(-N * u)
    logb = -N * u
    
    # Force default degeneracy
    gk_fn = lambda k: math.comb(d + k - 1, k)
    if logwk_fn is None:
        Ek_fn = lambda k: k
        logwk_fn = lambda k, logb: Ek_fn(k) * logb
    logZ = [-math.inf] * (n + 1)
    logZ[0] = 0.0
    def convolve_shell_logspace(logZ, g, logw, n):
        g = int(g)
        mmax = min(g, n)
        logC = [0.0] * (mmax + 1)
        logC[0] = 0.0
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

    shells_used = 0
    prev_logZn = -math.inf
    small_count = 0
    cumulative_capacity = 0
    for offset in range(safety_cap + 1):
        k = k_start + offset
        gk = gk_fn(k)
        logwk = logwk_fn(k, logb)
        convolve_shell_logspace(logZ, gk, logwk, n)
        shells_used += 1
        cumulative_capacity += int(gk)
        if cumulative_capacity >= n:
            curr = logZ[n]
            if prev_logZn != -math.inf and curr != -math.inf:
                delta_log = abs(curr - prev_logZn)
                if delta_log <= tol: small_count += 1
                else: small_count = 0
            prev_logZn = curr
            if small_count >= consecutive_small: break
    return logZ[n]

val1 = fermion_logZ_numeric(tau=2.0, N=200000, n=100, d=2, k_start=1)
val2 = fermion_logZ_numeric(tau=2.5, N=200000, n=100, d=2, k_start=1)
print((val1 - val2)/0.5)
