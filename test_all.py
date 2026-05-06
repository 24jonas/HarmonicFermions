import math
import numpy as np
import numba

@numba.njit
def logsumexp2(x, y):
    m = x if x > y else y
    return m + math.log(math.exp(x - m) + math.exp(y - m))

@numba.njit
def logsumexp_array(vals):
    m = np.max(vals)
    s = 0.0
    for i in range(len(vals)):
        s += math.exp(vals[i] - m)
    return m + math.log(s)

@numba.njit
def log1mexp_neg(x):
    return math.log(-math.expm1(x))

@numba.njit
def log_Q(n, b):
    logb = math.log(b)
    L_prev = np.zeros(n + 2, dtype=np.float64)
    L = np.zeros(n + 2, dtype=np.float64)

    for m in range(2, n + 2):
        L[m - 1] = logsumexp_array(L_prev[:m - 1])
        log1mbpow = log1mexp_neg((m - 1) * logb)   
        for i in range(m - 1, 0, -1):
            t1 = logb + L[i]                       
            t2 = (1 - i) * logb + log1mbpow + L_prev[i - 1]
            L[i - 1] = logsumexp2(t1, t2)
        temp = L_prev
        L_prev = L
        L = temp

    logH = L_prev[n]
    log_poch = 0.0
    for k in range(1, n + 1):
        log_poch += log1mexp_neg(k * logb)

    return (n * (n + 1) // 2) * logb - 2.0 * log_poch + logH

@numba.njit
def log_Q_all(n, b):
    if not (0.0 < b < 1.0):
        raise ValueError("b must lie in (0,1)")
    
    logb = math.log(b)
    L_prev = np.zeros(n + 2, dtype=np.float64)
    L = np.zeros(n + 2, dtype=np.float64)
    
    res = np.zeros(n + 1, dtype=np.float64)
    log_poch = 0.0
    
    for m in range(2, n + 2):
        L[m - 1] = logsumexp_array(L_prev[:m - 1])
        log1mbpow = log1mexp_neg((m - 1) * logb)   
        for i in range(m - 1, 0, -1):
            t1 = logb + L[i]                       
            t2 = (1 - i) * logb + log1mbpow + L_prev[i - 1]
            L[i - 1] = logsumexp2(t1, t2)
            
        temp = L_prev
        L_prev = L
        L = temp
        
        num = m - 1
        logH_num = L_prev[num]
        
        log_poch += log1mexp_neg(num * logb)
        res[num] = (num * (num + 1) / 2.0) * logb - 2.0 * log_poch + logH_num
        
    return res

n = 100
b = 0.5
res1 = log_Q(n, b)
res2 = log_Q_all(n, b)

print(f"res1: {res1}")
print(f"res2[-1]: {res2[-1]}")
print(f"Matches: {abs(res1 - res2[-1]) < 1e-10}")

res3 = log_Q(n-1, b)
print(f"res3: {res3}")
print(f"res2[-2]: {res2[-2]}")
print(f"Matches: {abs(res3 - res2[-2]) < 1e-10}")

res4 = log_Q(1, b)
print(f"res4: {res4}")
print(f"res2[1]: {res2[1]}")
print(f"Matches: {abs(res4 - res2[1]) < 1e-10}")
