particles = 0
energy = 0
for k in range(1, 30):
    to_add = min(k, 100 - particles)
    energy += to_add * k
    particles += to_add
    if particles == 100:
        break
print(f"Energy for 100 particles in 2D HO: {energy}")
