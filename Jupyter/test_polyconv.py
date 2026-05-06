import math
import numpy as np

def level_degeneracy_numba(d, k):
    if k == 0: return 1.0
    if d == 1: return 1.0
    res = 1.0
    for i in range(1, k + 1):
        res = res * (d + k - i) / i
    return res

print(level_degeneracy_numba(3, 2))
print(math.comb(3+2-1, 2))
