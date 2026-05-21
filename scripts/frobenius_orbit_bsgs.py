#!/usr/bin/env python3
"""
REFINED IDEA: Frobenius-Orbit BSGS

The key insight: c_{d+80}[i] = c_d[i] * p^i mod ell (VERIFIED).

This means the 96 diagonals group into 16 Frobenius orbits of size 6:
  Orbit 0: diagonals {0, 80, 64, 48, 32, 16}  (d, d+80, d+160%96, ...)
  Orbit 1: diagonals {1, 81, 65, 49, 33, 17}
  ...
  Orbit 15: diagonals {15, 95, 79, 63, 47, 31}

For each orbit, we need ONE rotation σ_0^d(x), and then the contributions
of all 6 diagonals in the orbit can be computed by multiplying by different
plaintext constants. NO additional key-switches needed!

In the standard BSGS for 96 diagonals: baby=10, giant=10 → 18 key-switches.
In the Frobenius-orbit BSGS for 16 effective diagonals: baby=4, giant=4 → 7 key-switches!

But wait - we need to check what the 16 orbit representatives' offsets are,
and whether they can be covered by a BSGS with fewer automorphisms.
"""
import numpy as np
import math

def power_mod(base, exp, mod):
    result = 1
    base = base % mod
    while exp > 0:
        if exp % 2 == 1:
            result = (result * base) % mod
        exp //= 2
        base = (base * base) % mod
    return result

ell = 97
p = 65537
D = 96
k_p = 80  # Frobenius shift: column j → column (j+80) mod 96

# Compute the 16 Frobenius orbits
orbits = []
visited = set()
for d in range(D):
    if d not in visited:
        orbit = []
        current = d
        for _ in range(6):  # orbit size = 6
            orbit.append(current)
            visited.add(current)
            current = (current + k_p) % D
        orbits.append(orbit)

print(f"Number of Frobenius orbits: {len(orbits)}")
print(f"Orbit size: 6 (= 96 / gcd(80, 96) = 96/16 = 6)")
print()

# The orbit representatives (smallest element of each orbit)
representatives = [min(orbit) for orbit in orbits]
representatives.sort()
print(f"Orbit representatives (16 values): {representatives}")
print()

# These are the ONLY rotation offsets we need!
# BSGS on these 16 offsets:
num_eff = len(representatives)
baby = int(math.ceil(math.sqrt(num_eff)))
giant = int(math.ceil(num_eff / baby))
print(f"Standard BSGS on {num_eff} offsets: baby={baby}, giant={giant} → {baby+giant-1} auts")
print(f"Current BSGS on 96 offsets: baby=10, giant=10 → 18 auts")
print(f"Reduction: {18 - (baby+giant-1)} fewer = {(18-(baby+giant-1))/18*100:.0f}%")
print()

# But the 16 representatives are NOT consecutive! They are: 0,1,2,...,15
# (since gcd(80,96)=16, the orbits are {d, d+16, d+32, d+48, d+64, d+80} for d=0..15)
# Wait, let me recheck:
print("Actual orbits:")
for i, orbit in enumerate(orbits[:5]):
    print(f"  Orbit {i}: {sorted(orbit)}")

# The representatives are 0,1,2,...,15 (consecutive!)
# BSGS on offsets {0,1,2,...,15}: baby=4, giant=4 → 7 auts
# Or baby=5, giant=4 → 8 auts (to cover 16 with 5*4=20 ≥ 16)

# Actually for 16 consecutive offsets: baby=4, giant=4 covers 4*4=16 exactly!
print(f"\nSince representatives are {{0,1,...,15}} (consecutive):")
print(f"  BSGS: baby=4 (σ^0,σ^1,σ^2,σ^3), giant=4 (σ^0,σ^4,σ^8,σ^12)")
print(f"  Total: 4+4-1 = 7 automorphisms!")
print()

# Cost model:
# Standard BSGS (D=96): 1 hoisting + 18 auts from precon + 96 MulAdds
# Frobenius-orbit BSGS: 1 hoisting + 7 auts from precon + 96 MulAdds
#   (still 96 MulAdds because we still multiply all 96 diagonals by their constants)
#
# Time estimate:
# Standard: 0.6 + 18*0.1 + 96*0.01 = 0.6 + 1.8 + 0.96 = 3.36s (per D=96 call)
# Wait, the actual time is 7.2s. Let me recalibrate.
#
# From profiling: mul_MatMul1DExec = 7.2s for D=96
# Components: BasicAutomorphPrecon(0.6s) + automorphs(151/6≈25 per call, 25*0.1=2.5s)
#   + MulAdd + giant-step key-switches
#
# Actually in HElib's BSGS:
#   - Baby steps: 10 hoisted automorphisms (from precon) = 10 * 0.1s = 1.0s
#   - Giant steps: 9 unhoisted automorphisms (full key-switch) = 9 * 0.5s = 4.5s
#   - MulAdd: 96 * 0.01s = 0.96s
#   - Precon: 0.6s
#   Total: 0.6 + 1.0 + 4.5 + 0.96 ≈ 7.06s ✓
#
# With Frobenius-orbit BSGS (7 auts total):
#   If we use baby=4, giant=4:
#   - Baby steps: 4 hoisted = 4 * 0.1s = 0.4s
#   - Giant steps: 3 unhoisted = 3 * 0.5s = 1.5s
#   - MulAdd: 96 * 0.01s = 0.96s (still 96 because each orbit member needs its own MulAdd)
#   - Precon: 0.6s
#   Total: 0.6 + 0.4 + 1.5 + 0.96 ≈ 3.46s
#   Speedup: 7.06 / 3.46 = 2.04x for the D=96 block!

print("=== COST MODEL ===")
print()
print("Current BSGS (D=96, baby=10, giant=10):")
print(f"  Precon: 0.6s")
print(f"  Baby steps (10 hoisted): 10 × 0.1s = 1.0s")
print(f"  Giant steps (9 unhoisted): 9 × 0.5s = 4.5s")
print(f"  MulAdd (96 diagonals): 96 × 0.01s = 0.96s")
print(f"  Total: 7.06s")
print()
print("Frobenius-Orbit BSGS (16 effective diags, baby=4, giant=4):")
print(f"  Precon: 0.6s")
print(f"  Baby steps (4 hoisted): 4 × 0.1s = 0.4s")
print(f"  Giant steps (3 unhoisted): 3 × 0.5s = 1.5s")
print(f"  MulAdd (96 diagonals, 6 per rotation): 96 × 0.01s = 0.96s")
print(f"  Total: 3.46s")
print(f"  Speedup: 7.06 / 3.46 = {7.06/3.46:.2f}x")
print()
print("Impact on full bootstrapping:")
print(f"  linear2 has 2 calls to D=96 MatMul (forward + inverse)")
print(f"  Current linear2: ~45s (includes D=29 dimension too)")
print(f"  D=96 portion: ~2 × 7.2s = 14.4s")
print(f"  With Frobenius-orbit: ~2 × 3.5s = 7.0s")
print(f"  Savings: ~7.4s on linear2")
print(f"  New linear2: ~45 - 7.4 = ~37.6s")
print(f"  Speedup on linear2: 45/37.6 = {45/37.6:.2f}x (16% faster)")
print()
print("=== IMPLEMENTATION APPROACH ===")
print()
print("Key change: modify MatMul1DExec to recognize Frobenius-orbit structure.")
print("Instead of iterating over 96 diagonals with BSGS(baby=10, giant=10),")
print("iterate over 16 orbit representatives with BSGS(baby=4, giant=4).")
print("For each representative d, compute 6 MulAdds using the same rotated ciphertext")
print("but with 6 different plaintext constants (one per orbit member).")
print()
print("The plaintext constants for orbit member d+k*80 are:")
print("  c_{d+k*80}[i] = c_d[i] * (p^i)^k mod ell")
print("  = c_d[i] * p^{ik} mod ell")
print()
print("This is a SINGLE-STAGE optimization (no multi-stage overhead)!")
print("It works within the existing BSGS framework.")
print("The only change is the diagonal grouping strategy.")
