#!/usr/bin/env python3
"""
Precise verification of the Good-Thomas M2 solving logic,
matching EXACTLY what the C++ code does.

The C++ code does:
  for d in range(D):
    k = ((d % N2) * inv_N1_mod_N2) % N2
    d1 = k * N1
    d2 = (d - d1 + D) % D
    if d2 % N2 != 0: continue
    for i in range(D):
      pos = (i + d2) % D
      exp = ((pos % N2) * (d1 / N1)) % N2
      inv_exp = (N2 - exp) % N2
      c1_inv = omega_N2^inv_exp
      c2_val = A_full[i][(i+d)%D] * c1_inv
      M2[i][(i+d2)%D] += c2_val

The issue: multiple values of d can map to the same d2!
When d1=0: d2=d (only if d%N2==0, i.e., d is multiple of 3)
When d1=32: d2=d-32 (only if (d-32)%N2==0)
When d1=64: d2=d-64 (only if (d-64)%N2==0)

For each d, exactly ONE (d1, d2) pair is selected.
But different d values CAN map to the same d2.
Example: d=0 → d1=0, d2=0; d=32 → d1=32, d2=0; d=64 → d1=64, d2=0
All three map to d2=0!

So M2[i][i%D] gets contributions from d=0, d=32, and d=64.
This is CORRECT - M2's diagonal at offset 0 should be the sum of
contributions from all d values that map to d2=0.

But wait - in the C++ code, M2[i][(i+d2)%D] += c2_val uses +=.
In NTL's Mat<RX>, += on RX should work correctly (polynomial addition).

Let me trace through a specific example to find the bug.
"""
import numpy as np

def power_mod(base, exp, mod):
    result = 1
    base = base % mod
    while exp > 0:
        if exp % 2 == 1:
            result = (result * base) % mod
        exp //= 2
        base = (base * base) % mod
    return result

def find_generator(m):
    phi = m - 1
    factors = []
    n = phi
    for p_f in range(2, int(n**0.5) + 1):
        if n % p_f == 0:
            factors.append(p_f)
            while n % p_f == 0:
                n //= p_f
    if n > 1:
        factors.append(n)
    for g in range(2, m):
        if all(power_mod(g, phi // f, m) != 1 for f in factors):
            return g
    return None

ell = 97
D = 96
N1 = 32
N2 = 3
g_ell = find_generator(ell)
cofactor = 523
cof_mod_ell = cofactor % ell

# Build A_full (scalar over F_ell)
A_full = np.zeros((D, D), dtype=np.int64)
for i in range(D):
    for j in range(D):
        point_j = (cof_mod_ell * power_mod(g_ell, j, ell)) % ell
        A_full[i][j] = power_mod(point_j, i, ell)

# omega_N2 = points[1]^N1
point_1 = (cof_mod_ell * power_mod(g_ell, 1, ell)) % ell
omega_N2 = power_mod(point_1, N1, ell)
print(f"omega_N2 = {omega_N2} (should have order 3)")
print(f"omega_N2^3 mod {ell} = {power_mod(omega_N2, 3, ell)} (should be 1)")

# Build M1
M1 = np.zeros((D, D), dtype=np.int64)
for d_idx in range(N2):
    d = d_idx * N1
    for i in range(D):
        exp = ((i % N2) * d_idx) % N2
        val = power_mod(omega_N2, exp, ell)
        M1[i][(i + d) % D] = val

# Build M2 using the EXACT same logic as C++
inv_N1_mod_N2 = pow(N1, -1, N2)  # inv(32, 3) = inv(2, 3) = 2
print(f"inv_N1_mod_N2 = {inv_N1_mod_N2}")

M2 = np.zeros((D, D), dtype=np.int64)
for d in range(D):
    k = ((d % N2) * inv_N1_mod_N2) % N2
    d1 = k * N1
    d2 = (d - d1 + D) % D
    if d2 % N2 != 0:
        print(f"SKIP: d={d}, k={k}, d1={d1}, d2={d2} (d2%N2={d2%N2})")
        continue
    for i in range(D):
        pos = (i + d2) % D
        exp = ((pos % N2) * (d1 // N1)) % N2
        inv_exp = (N2 - exp) % N2
        c1_inv = power_mod(omega_N2, inv_exp, ell)
        a_val = A_full[i][(i + d) % D]
        c2_val = (a_val * c1_inv) % ell
        M2[i][(i + d2) % D] = (M2[i][(i + d2) % D] + c2_val) % ell

# Verify M2 * M1 = A_full
product = np.zeros((D, D), dtype=np.int64)
for i in range(D):
    for j in range(D):
        s = 0
        for k in range(D):
            s = (s + M2[i][k] * M1[k][j]) % ell
        product[i][j] = s

errors = np.sum(product != A_full)
print(f"\nVerification M2*M1 = A_full: errors = {errors}")

if errors > 0:
    # Find first error
    for i in range(D):
        for j in range(D):
            if product[i][j] != A_full[i][j]:
                print(f"  First error at [{i}][{j}]: product={product[i][j]}, expected={A_full[i][j]}")
                break
        else:
            continue
        break
else:
    print("*** PYTHON LOGIC IS CORRECT ***")
    print("\nThe C++ bug must be in how polynomial arithmetic differs from scalar arithmetic.")
    print("Key difference: in C++, the entries are polynomials in F_p[X]/G(X),")
    print("not scalars in F_ell. The omega_N2 computation might differ.")
    print()
    print("In C++: omega_N2 = PowerMod(points[1], N1, G)")
    print("  points[1] = X^{cofactor * g^1} mod G = X^{523*5} mod G = X^2615 mod G")
    print("  omega_N2 = (X^2615)^32 mod G = X^{83680} mod G")
    print()
    print("The issue might be that (X^2615)^32 mod G is NOT a primitive 3rd root")
    print("in F_p[X]/G(X) for the INVERSE matrix case!")
    print()
    print("For the INVERSE matrix: A_full is inverted BEFORE the factorization.")
    print("The inverse of a Vandermonde is NOT necessarily a Vandermonde with")
    print("the same structure. So the Good-Thomas factorization that works for")
    print("the forward matrix may NOT work for the inverse!")
    print()
    print("THIS IS LIKELY THE BUG: we apply Good-Thomas to the INVERTED matrix,")
    print("but the factorization formula assumes Vandermonde structure.")
