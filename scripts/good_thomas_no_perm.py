#!/usr/bin/env python3
"""
Key insight: the permutation itself can be implemented as a MatMul1D
with exactly 1 non-zero entry per row (a permutation matrix).
A permutation matrix has at most D non-zero diagonals, but each diagonal
has only 1 non-zero position. In HElib, this is still expensive.

Better approach: DON'T fold permutations. Instead:
1. Implement Stage 1 as a SPARSE matrix (3 diagonals in Ruritanian space)
2. Implement Stage 2 as a SPARSE matrix (32 diagonals in CRT space)
3. The permutation between stages is handled by choosing the RIGHT
   diagonal offsets for each stage.

The real question is: in HElib's slot rotation framework, what rotation
offsets correspond to the "natural" operations of each stage?

In HElib, rotation by d means: output[i] = input[(i+d) mod D].
The generator for dim=0 has order 96.

Stage 1 (3-point DFT): operates on groups of 3 elements spaced 32 apart.
  - Needs rotations by 32 and 64 (to access the 3 elements in each group)
  - This gives 3 non-zero diagonals at offsets {0, 32, 64}

Stage 2 (32-point DFT): operates on groups of 32 elements spaced 3 apart.
  - Needs rotations by 3, 6, 9, ..., 93 (to access the 32 elements)
  - This gives 32 non-zero diagonals at offsets {0, 3, 6, ..., 93}

BUT: Stage 1 and Stage 2 operate in DIFFERENT index spaces!
Stage 1 input is in "natural" order, output is in "partially transformed" order.
Stage 2 input is Stage 1's output.

The question is: can we design the stages so that BOTH use the same
linear slot indexing (0..95)?

Answer: YES, if we absorb the permutation into the DIAGONAL CONSTANTS.

For Stage 1 with input in natural order:
  The 3-point DFT groups elements at positions {j, j+32, j+64} for j=0..31.
  output[i] = Σ_{d∈{0,32,64}} c1_d[i] * input[(i+d) mod 96]
  where c1_d[i] depends on which "group" position i is in.

For Stage 2 with input = Stage 1 output:
  The 32-point DFT groups elements at positions {j, j+3, j+6, ..., j+93} for j=0,1,2.
  output[i] = Σ_{d∈{0,3,6,...,93}} c2_d[i] * input[(i+d) mod 96]

The key question: does Stage2 * Stage1 = Step2Matrix?
This depends on whether the "natural" grouping is compatible.

Let me verify: in natural order (0..95), the 3-point groups are:
  Group j: {j, j+32, j+64} for j=0,...,31
  These are 32 groups of 3 elements each.

After the 3-point DFT, the output at position i encodes the k2-th
frequency component of the group containing i.

For Stage 2, the 32-point groups are:
  Group j: {j, j+3, j+6, ..., j+93} for j=0,1,2
  These are 3 groups of 32 elements each.

The issue: Stage 1 groups by "mod 32" (elements 32 apart),
Stage 2 groups by "mod 3" (elements 3 apart).
These are DIFFERENT groupings of the same 96 elements.

For the factorization to work WITHOUT an intermediate permutation,
we need the Stage 1 output to be "ready" for Stage 2's grouping.
This is exactly what Good-Thomas guarantees when gcd(32,3)=1!

Let me verify this directly by building the matrices with these
specific diagonal structures and checking if their product = Step2Matrix.
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


def main():
    ell = 97
    D = 96
    N1 = 32
    N2 = 3
    cofactor = 523

    g_ell = find_generator(ell)
    cof_mod_ell = cofactor % ell
    omega = g_ell

    # Build Step2Matrix
    A_full = np.zeros((D, D), dtype=np.int64)
    for i in range(D):
        for j in range(D):
            A_full[i][j] = power_mod(cof_mod_ell * power_mod(g_ell, j, ell), i, ell)

    # Approach: build Stage 1 and Stage 2 directly in linear indexing
    # using the "natural" diagonal offsets.
    #
    # Stage 1: 3-point DFT on groups {i, i+32, i+64}
    #   M1[i, j] ≠ 0 only if (j-i) mod 96 ∈ {0, 32, 64}
    #   M1[i, i] = c1_0[i]
    #   M1[i, (i+32)%96] = c1_32[i]  (but in MatMul1D convention: M[i,j] with j=(i+d)%D)
    #   M1[i, (i+64)%96] = c1_64[i]
    #
    # What should c1_d[i] be?
    # The 3-point DFT with root ω_3 = ω^32 = 35 (order 3 mod 97):
    #   V_3 = [[1, 1, 1], [1, ω_3, ω_3²], [1, ω_3², ω_3]]
    #
    # For element at position i, its "position within the 3-group" is:
    #   k2 = ? (which frequency component does output[i] represent?)
    #
    # In the Good-Thomas framework:
    #   Input index j has components: j mod 32 (N1-component), j mod 3 (N2-component)
    #   The 3-point DFT transforms the N2-component (j mod 3) to frequency k2.
    #
    # For the 3-point DFT applied to elements spaced 32 apart:
    #   The group containing position i consists of: i, i+32, i+64 (mod 96)
    #   These have N2-components (mod 3): i%3, (i+32)%3, (i+64)%3
    #   Since 32 ≡ 2 (mod 3) and 64 ≡ 1 (mod 3):
    #   N2-components: i%3, (i+2)%3, (i+1)%3
    #
    # So the 3 elements in the group have N2-values that are a PERMUTATION of {0,1,2}
    # (specifically: i%3, (i%3+2)%3, (i%3+1)%3)
    #
    # The 3-point DFT output at position i should be:
    #   output[i] = Σ_{n2=0}^{2} V_3[k2(i), n2] * input[element with N2=n2 in group]
    #
    # where k2(i) = i mod 3 (the output frequency is determined by position mod 3)
    #
    # The element with N2=n2 in the group of i is at position:
    #   i + 32*((n2 - i%3) mod 3)... let me think more carefully.

    # Actually, let me just try the direct approach:
    # Stage 1 diagonal d=0: c1_0[i] = 1 (identity component of 3-point DFT)
    # Stage 1 diagonal d=32: c1_32[i] = ω_3^{(i mod 3) * 1} ... ?
    #
    # This is getting complicated. Let me just try ALL possible 3-diagonal
    # matrices and see which one, combined with a 32-diagonal Stage 2,
    # gives the Step2Matrix.

    # Actually, the cleanest approach: compute the CORRECT Stage 1 and Stage 2
    # by solving for them.
    #
    # We know: Step2Matrix = M2 * M1
    # M1 has diagonals at {0, 32, 64}
    # M2 has diagonals at {0, 3, 6, ..., 93}
    #
    # M1[i,j] = c1_{(j-i)%96}[i] if (j-i)%96 ∈ {0,32,64}, else 0
    # M2[i,j] = c2_{(j-i)%96}[i] if (j-i)%96 ∈ {0,3,...,93}, else 0
    #
    # (M2*M1)[i,k] = Σ_j M2[i,j] * M1[j,k]
    #              = Σ_{d2 ∈ {0,3,...,93}} Σ_{d1 ∈ {0,32,64}}
    #                  c2_d2[i] * c1_d1[(i+d2)%96]  if k = (i+d2+d1)%96
    #
    # This means: A[i, (i+d)%96] = Σ over (d2,d1) with d2+d1≡d (mod 96)
    #                                c2_d2[i] * c1_d1[(i+d2)%96]
    #
    # For each diagonal d of A, the constraint is:
    #   A_d[i] = Σ_{d2: d2 mod 3 = d mod 3, d2 ∈ {0,3,...,93}}
    #              c2_d2[i] * c1_{(d-d2)%96}[(i+d2)%96]
    #   where (d-d2)%96 must be in {0, 32, 64}
    #
    # (d-d2) mod 96 ∈ {0,32,64} means d2 ∈ {d, d-32, d-64} mod 96
    # AND d2 must be in {0,3,6,...,93} (multiples of 3)
    #
    # So for each d, the contributing (d2, d1) pairs are:
    #   d2 = d mod 96 (if d%3==0), d1 = 0
    #   d2 = (d-32) mod 96 (if (d-32)%3==0), d1 = 32
    #   d2 = (d-64) mod 96 (if (d-64)%3==0), d1 = 64
    #
    # Since d2 must be a multiple of 3:
    #   d%3==0 → d2=d, d1=0 contributes
    #   (d-32)%3==0 → d≡32≡2 (mod 3) → d2=d-32, d1=32 contributes
    #   (d-64)%3==0 → d≡64≡1 (mod 3) → d2=d-64, d1=64 contributes
    #
    # For any d, exactly ONE of {d, d-32, d-64} is ≡0 mod 3
    # (since 32≡2, 64≡1 mod 3, so {d, d-2, d-1} mod 3 covers all residues)
    #
    # Wait: d%3, (d-32)%3 = (d-2)%3 = (d+1)%3, (d-64)%3 = (d-1)%3 = (d+2)%3
    # So {d%3, (d+1)%3, (d+2)%3} = {0,1,2}. Exactly one is 0.
    #
    # This means each diagonal d of A_full is produced by exactly ONE (d2, d1) pair!
    # This is the key property of Good-Thomas (no twiddle factors needed).
    #
    # Therefore:
    #   If d%3 == 0: A_d[i] = c2_d[i] * c1_0[(i+d)%96]
    #   If d%3 == 1: (d-64)%3=0, so d2=(d-64)%96, d1=64
    #                A_d[i] = c2_{(d-64)%96}[i] * c1_64[(i+(d-64)%96)%96]
    #   If d%3 == 2: (d-32)%3=0, so d2=(d-32)%96, d1=32
    #                A_d[i] = c2_{(d-32)%96}[i] * c1_32[(i+(d-32)%96)%96]

    # This gives us a system we can solve!
    # If we FIX c1 (the 3-point DFT constants), we can solve for c2.
    #
    # Natural choice for c1 (3-point DFT):
    #   c1_0[i] = 1 for all i (sum of all 3 inputs)
    #   c1_32[i] = ω_3^{f(i)} for some function f
    #   c1_64[i] = ω_3^{2*f(i)} for some function f
    #
    # The function f(i) determines which "frequency" output[i] represents.
    # Natural choice: f(i) = i mod 3 (output position mod 3 = frequency index)

    omega_3 = power_mod(omega, N1, ell)  # ω^32 = 35, order 3
    omega_32 = power_mod(omega, N2, ell)  # ω^3 = 28, order 32

    print(f"ω_3 = {omega_3} (order 3), ω_32 = {omega_32} (order 32)")

    # Try: c1_0[i] = 1, c1_32[i] = ω_3^{i%3}, c1_64[i] = ω_3^{2*(i%3)}
    c1 = {0: np.ones(D, dtype=np.int64),
           32: np.array([power_mod(omega_3, i % N2, ell) for i in range(D)], dtype=np.int64),
           64: np.array([power_mod(omega_3, 2 * (i % N2), ell) for i in range(D)], dtype=np.int64)}

    # Build M1
    M1 = np.zeros((D, D), dtype=np.int64)
    for d in [0, 32, 64]:
        for i in range(D):
            j = (i + d) % D
            M1[i][j] = c1[d][i]

    # Now solve for c2 from: A_d[i] = c2_{d2}[i] * c1_{d1}[(i+d2)%96]
    c2 = {}
    for d2 in range(0, D, N2):  # d2 = 0, 3, 6, ..., 93
        c2[d2] = np.zeros(D, dtype=np.int64)

    for d in range(D):
        if d % 3 == 0:
            d2, d1 = d, 0
        elif d % 3 == 1:
            d2, d1 = (d - 64) % D, 64
        else:  # d % 3 == 2
            d2, d1 = (d - 32) % D, 32

        for i in range(D):
            A_d_i = A_full[i][(i + d) % D]
            c1_val = c1[d1][(i + d2) % D]
            if c1_val == 0:
                if A_d_i != 0:
                    print(f"ERROR: A_d[{i}] != 0 but c1 = 0 at d={d}")
                continue
            c1_inv = power_mod(int(c1_val), ell - 2, ell)
            c2[d2][i] = (A_d_i * c1_inv) % ell

    # Build M2
    M2 = np.zeros((D, D), dtype=np.int64)
    for d2 in range(0, D, N2):
        for i in range(D):
            j = (i + d2) % D
            M2[i][j] = c2[d2][i]

    # Verify M2 * M1 = A_full
    product = M2 @ M1 % ell
    errors = np.sum(product != A_full)
    print(f"\nVerification: M2 * M1 = Step2Matrix? Errors = {errors}")

    if errors == 0:
        print("\n*** SUCCESS: 2-stage decomposition with NO intermediate permutation! ***")
        print(f"\nStage 1: 3 non-zero diagonals at offsets {{0, 32, 64}}")
        print(f"  Automorphisms: 2 (σ^32, σ^64)")
        print(f"\nStage 2: 32 non-zero diagonals at offsets {{0, 3, 6, ..., 93}}")
        print(f"  Automorphisms with BSGS: baby=6(σ^3), giant=6(σ^18) → 11")
        print(f"\nTotal: 2 + 11 = 13 automorphisms (vs 18 for current BSGS)")
        print(f"Reduction: 27.8%")
        print(f"\nNoise: only 2 sequential key-switch stages (vs 13 for Rader97)")
        print(f"\n=== Stage 1 diagonal constants (for C++ implementation) ===")
        print(f"c1_0[i] = 1 for all i")
        print(f"c1_32[i] = ω_3^(i mod 3) = {omega_3}^(i mod 3)")
        print(f"c1_64[i] = ω_3^(2*(i mod 3)) = {omega_3}^(2*(i mod 3))")
        print(f"\nω_3 values: ω_3^0=1, ω_3^1={omega_3}, ω_3^2={power_mod(omega_3,2,ell)}")

        # Verify c2 structure
        print(f"\n=== Stage 2 diagonal constants ===")
        print(f"32 diagonals at offsets 0, 3, 6, ..., 93")
        # Check if c2 has any pattern
        print(f"Sample c2_0: {c2[0][:6]}...")
        print(f"Sample c2_3: {c2[3][:6]}...")
        print(f"Sample c2_6: {c2[6][:6]}...")
    else:
        print(f"FAILED: {errors} mismatches. Trying different c1 choices...")
        # Try f(i) = (i * inv_N1_mod_N2) % N2 or other mappings
        # The issue might be the twiddle factor assignment


if __name__ == "__main__":
    main()
