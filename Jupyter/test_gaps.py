import math

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
        
max_gap = 0
current_gap = 0
for k in range(1, 2000):
    g = degeneracy_box(2, k)
    if g == 0:
        current_gap += 1
        if current_gap > max_gap:
            max_gap = current_gap
    else:
        if current_gap >= 8:
            print(f"Gap of {current_gap} ends at {k}")
        current_gap = 0
