
def energy(n, w, verbose=False):
    """
    Computes the exact energy for n fermions in a 2D harmonic trap, 
    scaled by the effective relative frequency w.
    """
    # 1. Compute the standard base energy (w=1) using your exact state-filling logic
    cur_lvl = 0
    total_eng = 0
    m = n
    while n > 0:
        for i in range(0, cur_lvl):
            if n > 0:
                total_eng += cur_lvl - 1
                n -= 1
            else:
                break
        cur_lvl += 1
        
    base_energy = total_eng + m   
    
    if verbose:
        print("Highest shell: ", cur_lvl)
        
    # 2. Scale the relative energy component by w
    # For a 2D harmonic oscillator, the 1-particle CM energy is exactly 1.0
    e_1 = 1.0 
    
    true_energy_w = e_1 + w * (base_energy - e_1)
    
    return true_energy_w



print(energy(10000, 0.5))