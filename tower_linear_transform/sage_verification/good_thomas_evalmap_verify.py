#!/usr/bin/env python3
"""
Verify Good-Thomas PFA factorization of the D=96 Step2Matrix.

Key idea: 96 = 32 × 3, gcd(32,3) = 1.
Under CRT re-indexing, the Vandermonde matrix V[i,j] = ω_j^i
factors as V_32 ⊗ V_3 (tensor product).

This means the 96-point linear transform can be computed as:
  Stage 1: 3-point transform (2 automorphisms, 3 non-zero diagonals)
  Stage 2: 32-point transform (10 automorphisms via BSGS)
  Total: 12 automorphisms vs current 18 (33% reduction)

Crucially: only 2 sequential stages → manageable noise (vs Rader97's 13 stages).
"""

import numpy as np
from itertools import product


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
    """Find a generator of (Z/mZ)* for prime m."""
    phi = m - 1
    # Factor phi
    factors = []
    n = phi
    for p in range(2, int(n**0.5) + 1):
        if n % p == 0:
            factors.append(p)
            while n % p == 0:
                n //= p
    if n > 1:
        factors.append(n)

    for g in range(2, m):
        is_gen = True
        for f in factors:
            if power_mod(g, phi // f, m) == 1:
                is_gen = False
                break
        if is_gen:
            return g
    return None


def build_step2_matrix_Fp(p, ell, D):
    """
    Build the Step2Matrix over F_p (scalar version, ignoring polynomial ring).

    For our parameter set:
      p = 65537, m = 50731 = 97 × 523
      For dim=0 (factor 97): D = ord_97(p) = 96, cofactor = 523

    The matrix is: A[i,j] = (cofactor * reps[j])^i mod ell
    where reps[j] = g^j mod ell, g = generator of (Z/ell Z)*.

    But since we work mod p (in the slot polynomial), the actual computation is:
    A[i,j] = point_j^i mod p
    where point_j = (cofactor * reps[j]) mod ell, lifted to F_p.

    For the Vandermonde structure analysis, we just need:
    A[i,j] = omega_j^i mod p
    where omega_j = g^j mod ell (the evaluation points).
    """
    # Find generator of (Z/ell Z)*
    g_ell = find_generator(ell)

    # In HElib: g = InvMod(ZmStarGen(dim) % m, m)
    # For simplicity, use g_ell directly as the generator
    # reps[j] = g^j mod ell
    cofactor = 50731 // ell  # = 523 for ell=97

    # Build evaluation points
    # point_j = (cofactor * g^j mod ell) mod p
    # Since ell=97 is prime and cofactor=523, cofactor mod ell = 523 mod 97 = 523 - 5*97 = 523-485 = 38
    cof_mod_ell = cofactor % ell

    points = np.zeros(D, dtype=np.int64)
    for j in range(D):
        rep_j = power_mod(g_ell, j, ell)
        points[j] = (cof_mod_ell * rep_j) % ell

    # Build Vandermonde matrix mod p
    A = np.zeros((D, D), dtype=np.int64)
    for j in range(D):
        A[0][j] = 1
        for i in range(1, D):
            A[i][j] = (A[i-1][j] * points[j]) % p

    return A, points


def crt_index(a, b, D1=32, D2=3):
    """CRT: map (a mod D1, b mod D2) -> index mod D1*D2."""
    # Find x such that x ≡ a (mod D1) and x ≡ b (mod D2)
    # x = a * D2 * inv(D2, D1) + b * D1 * inv(D1, D2)
    D = D1 * D2
    inv_D2_mod_D1 = power_mod(D2, D1 - 2, D1) if D1 > 2 else 1  # D1=32 is power of 2, need different
    # For D1=32, D2=3: inv(3, 32) = 11 (since 3*11=33≡1 mod 32)
    inv_D2_mod_D1 = pow(D2, -1, D1)
    inv_D1_mod_D2 = pow(D1, -1, D2)
    x = (a * D2 * inv_D2_mod_D1 + b * D1 * inv_D1_mod_D2) % D
    return x


def build_crt_permutation(D1=32, D2=3):
    """Build CRT permutation: perm[linear_idx] = crt_idx."""
    D = D1 * D2
    perm = np.zeros(D, dtype=int)
    inv_perm = np.zeros(D, dtype=int)
    for a in range(D1):
        for b in range(D2):
            linear = a * D2 + b  # row-major (a, b)
            crt = crt_index(a, b, D1, D2)
            perm[linear] = crt
            inv_perm[crt] = linear
    return perm, inv_perm


def check_tensor_factorization(A, p, D1=32, D2=3):
    """
    Check if A (under CRT re-indexing) factors as V_D1 ⊗ V_D2.

    Under CRT indexing, if A is a pure Vandermonde V[i,j] = ω^{ij}:
      V[CRT(c,d), CRT(a,b)] = ω^{CRT(c,d) * CRT(a,b)}

    For Good-Thomas to work, we need:
      V = P^{-1} (V_D1 ⊗ V_D2) P
    where P is the CRT permutation.
    """
    D = D1 * D2
    perm, inv_perm = build_crt_permutation(D1, D2)

    # Reindex A: B[linear_i, linear_j] = A[perm[linear_i], perm[linear_j]]
    B = np.zeros((D, D), dtype=np.int64)
    for i in range(D):
        for j in range(D):
            B[i][j] = A[perm[i]][perm[j]]

    # Check if B factors as tensor product
    # B[(a,b), (c,d)] should equal V1[a,c] * V2[b,d] mod p
    # where (a,b) is the row-major index: linear_i = a*D2 + b

    # Extract V1 and V2 candidates
    # V1[a,c] = B[(a,0), (c,0)] (fix b=d=0)
    V1 = np.zeros((D1, D1), dtype=np.int64)
    for a in range(D1):
        for c in range(D1):
            V1[a][c] = B[a * D2 + 0][c * D2 + 0]

    # V2[b,d] = B[(0,b), (0,d)] (fix a=c=0)
    V2 = np.zeros((D2, D2), dtype=np.int64)
    for b in range(D2):
        for d in range(D2):
            V2[b][d] = B[0 * D2 + b][0 * D2 + d]

    # Verify: B[(a,b), (c,d)] == V1[a,c] * V2[b,d] mod p
    errors = 0
    for a in range(D1):
        for b in range(D2):
            for c in range(D1):
                for d in range(D2):
                    expected = (V1[a][c] * V2[b][d]) % p
                    actual = B[a * D2 + b][c * D2 + d]
                    if expected != actual:
                        errors += 1

    return errors, V1, V2


def analyze_diagonal_structure(A, D, p):
    """Analyze which diagonals are non-zero and their relationships."""
    diags = {}
    for d in range(D):
        diag = np.array([(A[i][(i + d) % D]) for i in range(D)], dtype=np.int64)
        if np.any(diag != 0):
            diags[d] = diag

    print(f"Total non-zero diagonals: {len(diags)} / {D}")
    return diags


def analyze_sparse_stage_structure(V1, V2, D1, D2, p):
    """
    Analyze the two-stage structure:
    Stage 1 (3-point): operates on diagonals {0, 32, 64}
    Stage 2 (32-point): operates on diagonals {0, 3, 6, ..., 93}
    """
    print(f"\n=== Stage 1: {D2}-point transform ===")
    print(f"Non-zero diagonals of V2 ({D2}x{D2}):")
    for d in range(D2):
        diag = [(V2[i][(i + d) % D2]) for i in range(D2)]
        print(f"  diagonal {d} (offset {d*D1} in full matrix): {diag}")

    print(f"\n=== Stage 2: {D1}-point transform ===")
    # Count non-zero diagonals of V1
    nz_diags_V1 = 0
    for d in range(D1):
        diag = [(V1[i][(i + d) % D1]) for i in range(D1)]
        if any(x != 0 for x in diag):
            nz_diags_V1 += 1
    print(f"Non-zero diagonals of V1 ({D1}x{D1}): {nz_diags_V1}")

    # BSGS estimate for V1
    import math
    baby = int(math.ceil(math.sqrt(nz_diags_V1)))
    giant = int(math.ceil(nz_diags_V1 / baby))
    bsgs_auts = baby + giant - 1
    print(f"BSGS for {D1}-point: baby={baby}, giant={giant}, auts={bsgs_auts}")

    # Stage 1 automorphisms
    stage1_auts = D2 - 1  # 3-point needs 2 automorphisms
    print(f"\nStage 1 automorphisms: {stage1_auts}")
    print(f"Stage 2 automorphisms: {bsgs_auts}")
    print(f"Total: {stage1_auts + bsgs_auts}")
    print(f"Current BSGS for D={D1*D2}: baby=10, giant=10, auts=18")
    print(f"Reduction: {18 - (stage1_auts + bsgs_auts)} fewer automorphisms ({(18-(stage1_auts+bsgs_auts))/18*100:.1f}%)")


def main():
    p = 65537
    ell = 97  # prime factor of m=50731
    D = 96    # ord_97(65537) = 96
    D1 = 32
    D2 = 3

    print(f"Parameters: p={p}, ell={ell}, D={D}={D1}×{D2}")
    print(f"gcd({D1},{D2}) = {np.gcd(D1, D2)}")
    print(f"{D1} | p-1 = {p-1}? {(p-1) % D1 == 0}")
    print(f"{D2} | p-1 = {p-1}? {(p-1) % D2 == 0}")
    print()

    # Build the Step2Matrix (scalar version over F_p)
    print("Building Step2Matrix (Vandermonde over F_p)...")
    A, points = build_step2_matrix_Fp(p, ell, D)
    print(f"Matrix size: {D}×{D}")
    print(f"First few evaluation points: {points[:5]}")

    # Analyze diagonal structure
    print("\n=== Diagonal Analysis ===")
    diags = analyze_diagonal_structure(A, D, p)

    # Check tensor factorization
    print("\n=== Good-Thomas Tensor Factorization Check ===")
    errors, V1, V2 = check_tensor_factorization(A, p, D1, D2)
    print(f"Factorization errors: {errors} / {D*D}")

    if errors == 0:
        print("*** PERFECT TENSOR FACTORIZATION! ***")
        print("Step2Matrix = P^{-1} (V_32 ⊗ V_3) P under CRT re-indexing")
        analyze_sparse_stage_structure(V1, V2, D1, D2, p)
    else:
        print(f"Tensor factorization has {errors} mismatches.")
        print("Checking if it's a 'twisted' tensor product (with twiddle factors)...")

        # Check if B[(a,b),(c,d)] = V1[a,c] * V2[b,d] * T[a,b,c,d] for some twiddle T
        perm, inv_perm = build_crt_permutation(D1, D2)
        B = np.zeros((D, D), dtype=np.int64)
        for i in range(D):
            for j in range(D):
                B[i][j] = A[perm[i]][perm[j]]

        # Check row-wise: for fixed (a,b), is B[(a,b), :] a scaled version of V1[a,:] ⊗ V2[b,:]?
        consistent_rows = 0
        for a in range(min(D1, 4)):
            for b in range(D2):
                row = B[a * D2 + b]
                # Expected: V1[a,c] * V2[b,d] * twiddle(a,b)
                ref_row = np.array([(V1[a][c] * V2[b][d]) % p
                                    for c in range(D1) for d in range(D2)], dtype=np.int64)
                # Check if row is a scalar multiple of ref_row
                if ref_row[0] != 0 and row[0] != 0:
                    ratio = (row[0] * pow(int(ref_row[0]), -1, p)) % p
                    scaled = (ref_row * ratio) % p
                    if np.array_equal(row, scaled):
                        consistent_rows += 1

        print(f"Rows consistent with scalar-twisted tensor: {consistent_rows}/{min(D1,4)*D2}")


if __name__ == "__main__":
    main()
