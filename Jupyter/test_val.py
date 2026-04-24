import math

def logaddexp(a, b):
    if a == -math.inf:
        return b
    if b == -math.inf:
        return a
    if a > b:
        return a + math.log1p(math.exp(b - a))
    else:
        return b + math.log1p(math.exp(a - b))

def degeneracy_box(d, k):
    if d == 1:
        root = int(math.sqrt(k))
        if root * root == k and root > 0:
            return 1
        return 0
    elif d == 2:
        count = 0
        max_n = int(math.sqrt(k))
        for nx in range(1, max_n + 1):
            ny_sq = k - nx**2
            if ny_sq > 0:
                ny = int(math.sqrt(ny_sq))
                if ny**2 == ny_sq:
                    count += 1
        return count
    else:
        count = 0
        max_n = int(math.sqrt(k))
        for nx in range(1, max_n + 1):
            count += degeneracy_box(d - 1, k - nx**2)
        return count

def log_binom(n, k):
    if k < 0 or k > n:
        return -math.inf
    return math.lgamma(n + 1) - math.lgamma(k + 1) - math.lgamma(n - k + 1)

def fermion_logZ_numeric(
    tau,
    N,
    n,
    d=None,
    max_shell=None,
    tol=1e-4,
    consecutive_small=8,
    safety_cap=100000,
    return_all=False,
    gk_fn=None,
    Ek_fn=None,
    logwk_fn=None,
    k_start=0,
):
    if N <= 0: raise ValueError("N must be a positive integer.")
    if n < 0: raise ValueError("n must be a nonnegative integer.")
    if tau < 0: raise ValueError("tau must be nonnegative.")
    if max_shell is not None and max_shell < k_start: raise ValueError("max_shell must be >= k_start.")
    if k_start < 0: raise ValueError("k_start must be nonnegative.")

    if gk_fn is None and d is None: raise ValueError("Provide d for the default degeneracy, or provide gk_fn.")
    if d is not None and d <= 0: raise ValueError("d must be a positive integer.")

    epsilon = tau / N
    zeta = 1.0 + 0.5 * epsilon * epsilon
    u = math.log(zeta + math.sqrt(zeta * zeta - 1.0))
    b = math.exp(-N * u)
    logb = -N * u  # log(b)

    if gk_fn is None:
        gk_fn = lambda k: degeneracy_box(d, k)

    if logwk_fn is None:
        if Ek_fn is None:
            Ek_fn = lambda k: k
        logwk_fn = lambda k, logb: Ek_fn(k) * logb

    logZ = [-math.inf] * (n + 1)
    logZ[0] = 0.0

    if n == 0:
        info = { "epsilon": epsilon, "zeta": zeta, "u": u, "b": b, "logb": logb, "shells_used": 0, "k_start": k_start, }
        return (0.0, logZ, info) if return_all else 0.0

    def convolve_shell_logspace(logZ, g, logw, n):
        if g < 0 or int(g) != g: raise ValueError(f"g_k must be a nonnegative integer, got {g}.")
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

    if max_shell is not None:
        shells_used = 0
        for k in range(k_start, max_shell + 1):
            gk = gk_fn(k)
            logwk = logwk_fn(k, logb)
            convolve_shell_logspace(logZ, gk, logwk, n)
            shells_used += 1
    else:
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
                    if delta_log <= tol:
                        small_count += 1
                    else:
                        small_count = 0
                prev_logZn = curr
                if small_count >= consecutive_small:
                    break
        else:
            raise RuntimeError("Adaptive truncation did not converge.")

    if return_all:
        info = { "epsilon": epsilon, "zeta": zeta, "u": u, "b": b, "logb": logb, "shells_used": shells_used, "k_start": k_start, }
        return logZ[n], logZ, info
    return logZ[n]

n = 100
val1 = fermion_logZ_numeric(tau=2.0, N=200000, n=n, d=2, k_start=1)
val2 = fermion_logZ_numeric(tau=2.5, N=200000, n=n, d=2, k_start=1)

print("logZ(2.0) =", val1)
print("logZ(2.5) =", val2)
print("val =", (val1 - val2) / 0.5)

