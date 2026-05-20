#!/usr/bin/env python3
"""
Correct Good-Thomas PFA verification for D=96 = 32 × 3.

The Good-Thomas PFA works by RE-INDEXING rows and columns via CRT,
then the DFT matrix becomes a tensor product.

Specifically: if we define the CRT bijection
  φ: Z/96Z → Z/32Z × Z/3Z,  φ(n) = (n mod 32, n mod 3)

Then the DFT matrix V[i,j] = ω^{ij} satisfies:
  V[φ^{-1}(a,b), φ^{-1}(c,d)] = ω^{φ^{-1}(a,b) * φ^{-1}(c,d)}

The key property of Good-Thomas is that under the INPUT permutation
σ_in(j) = j mod 32 * 3 + j mod 3 * 32 (Ruritanian map) and
OUTPUT permutation, the DFT becomes a tensor product WITHOUT twiddle factors.

Let me implement this correctly.
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
    p = 65537
    ell = 97
    D = 96
    N1 = 32
    N2 = 3

    g_ell = find_generator(ell)
    print(f"Generator of (Z/{ell}Z)*: g = {g_ell}, order = {D}")

    # Build DFT matrix: V[i,j] = g^{ij} mod ell
    omega = g_ell  # primitive 96th root of unity mod 97
    V = np.zeros((D, D), dtype=np.int64)
    for i in range(D):
        for j in range(D):
            V[i][j] = power_mod(omega, (i * j) % D, ell)

    # Good-Thomas PFA:
    # Input permutation (column): σ_in(j) maps (j1, j2) -> j
    #   where j1 = j mod N1, j2 = j mod N2
    #   j = N2 * j1 + N1 * j2  (mod N1*N2)  [Ruritanian map, no mod needed since gcd=1]
    #   Actually: j = (N2 * inv(N2,N1)) * j1 * N1_part + ...
    #
    # Standard Good-Thomas:
    #   Input map:  n = n1*N2 + n2*N1 mod N  (where n1 ∈ Z/N1, n2 ∈ Z/N2)
    #   Output map: k = k1*N2*inv(N2,N1) + k2*N1*inv(N1,N2) mod N  (CRT reconstruction)
    #
    # Then: V[k, n] = ω^{kn} and under these maps:
    #   kn mod N = [k1*N2*inv(N2,N1) + k2*N1*inv(N1,N2)] * [n1*N2 + n2*N1] mod N
    #
    # Expanding:
    #   = k1*N2*inv(N2,N1)*n1*N2 + k1*N2*inv(N2,N1)*n2*N1
    #     + k2*N1*inv(N1,N2)*n1*N2 + k2*N1*inv(N1,N2)*n2*N1  mod N
    #
    # Since N = N1*N2:
    #   k1*N2²*inv(N2,N1)*n1 mod N: N2²*inv(N2,N1) mod N1 = N2*1 = N2...
    #   This is getting complex. Let me just verify numerically.

    # Standard Good-Thomas input permutation: n -> (n mod N1, n mod N2)
    # Inverse: (n1, n2) -> n1 * N2 * pow(N2, -1, N1) + n2 * N1 * pow(N1, -1, N2) mod N
    inv_N2_N1 = pow(N2, -1, N1)  # inv(3, 32) = 11
    inv_N1_N2 = pow(N1, -1, N2)  # inv(32, 3) = 2

    def crt_reconstruct(n1, n2):
        """(n1 mod N1, n2 mod N2) -> n mod N"""
        return (n1 * N2 * inv_N2_N1 + n2 * N1 * inv_N1_N2) % D

    def ruritanian(n1, n2):
        """Good-Thomas input map: (n1, n2) -> n1*N2 + n2*N1 mod N"""
        return (n1 * N2 + n2 * N1) % D

    # Verify both are bijections
    crt_vals = set()
    rur_vals = set()
    for n1 in range(N1):
        for n2 in range(N2):
            crt_vals.add(crt_reconstruct(n1, n2))
            rur_vals.add(ruritanian(n1, n2))
    assert len(crt_vals) == D, "CRT not bijective!"
    assert len(rur_vals) == D, "Ruritanian not bijective!"
    print(f"Both CRT and Ruritanian maps are bijections on Z/{D}Z ✓")

    # Good-Thomas theorem:
    # If we permute COLUMNS by Ruritanian and ROWS by CRT:
    #   P_out * V * P_in = V_N1 ⊗ V_N2
    # where:
    #   P_in permutes column j to position (j mod N1, j mod N2) via Ruritanian^{-1}
    #   P_out permutes row i to position (i mod N1, i mod N2) via CRT^{-1}
    #
    # More precisely:
    #   V_reindexed[(k1,k2), (n1,n2)] = V[crt(k1,k2), rur(n1,n2)]
    #                                   = ω^{crt(k1,k2) * rur(n1,n2)}
    #
    # We need: crt(k1,k2) * rur(n1,n2) mod N = ?
    # crt(k1,k2) = k1*N2*inv(N2,N1) + k2*N1*inv(N1,N2)
    # rur(n1,n2) = n1*N2 + n2*N1
    #
    # Product mod N:
    # = [k1*N2*inv(N2,N1) + k2*N1*inv(N1,N2)] * [n1*N2 + n2*N1] mod N1*N2
    #
    # Term 1: k1*N2*inv(N2,N1) * n1*N2 = k1*n1*N2²*inv(N2,N1)
    #   mod N1: N2²*inv(N2,N1) mod N1 = N2*(N2*inv(N2,N1) mod N1) = N2*1 = N2...
    #   Actually mod N=N1*N2: this equals k1*n1*N2²*inv(N2,N1) mod N1*N2
    #   Since N2*inv(N2,N1) ≡ 1 mod N1, we have N2*inv(N2,N1) = 1 + t*N1 for some t
    #   So N2²*inv(N2,N1) = N2*(1+t*N1) = N2 + t*N1*N2 = N2 + t*N ≡ N2 mod N
    #   Wait no: N2*inv(N2,N1) ≡ 1 (mod N1), so N2*inv(N2,N1) = 1 + s*N1
    #   Then k1*n1*N2²*inv(N2,N1) = k1*n1*N2*(1+s*N1) = k1*n1*N2 + k1*n1*s*N mod N
    #   ≡ k1*n1*N2 mod N
    #
    # Term 2: k1*N2*inv(N2,N1) * n2*N1 = k1*n2*N1*N2*inv(N2,N1) = k1*n2*N*inv(N2,N1) ≡ 0 mod N
    #
    # Term 3: k2*N1*inv(N1,N2) * n1*N2 = k2*n1*N1*N2*inv(N1,N2) = k2*n1*N*inv(N1,N2) ≡ 0 mod N
    #
    # Term 4: k2*N1*inv(N1,N2) * n2*N1 = k2*n2*N1²*inv(N1,N2)
    #   Similarly: N1*inv(N1,N2) ≡ 1 mod N2, so N1²*inv(N1,N2) ≡ N1 mod N
    #   So this term ≡ k2*n2*N1 mod N
    #
    # Total: crt(k1,k2) * rur(n1,n2) ≡ k1*n1*N2 + k2*n2*N1 mod N
    #
    # Therefore: ω^{crt(k1,k2)*rur(n1,n2)} = ω^{k1*n1*N2 + k2*n2*N1}
    #          = ω^{k1*n1*N2} * ω^{k2*n2*N1}
    #          = (ω^{N2})^{k1*n1} * (ω^{N1})^{k2*n2}
    #
    # ω^{N2} = ω^3 has order N1=32 (primitive 32nd root)
    # ω^{N1} = ω^32 has order N2=3 (primitive 3rd root)
    #
    # So: V_reindexed[(k1,k2),(n1,n2)] = (ω^3)^{k1*n1} * (ω^32)^{k2*n2}
    #                                    = V_32[k1,n1] * V_3[k2,n2]
    #
    # This IS the tensor product V_32 ⊗ V_3 !!!

    print("\n=== Verifying Good-Thomas PFA ===")
    print(f"ω = g = {omega} (primitive 96th root mod {ell})")
    omega_N2 = power_mod(omega, N2, ell)  # ω^3, order 32
    omega_N1 = power_mod(omega, N1, ell)  # ω^32, order 3
    print(f"ω^{N2} = {omega_N2}, order = {96//np.gcd(N2, 96)} (should be {N1})")
    print(f"ω^{N1} = {omega_N1}, order = {96//np.gcd(N1, 96)} (should be {N2})")

    # Verify
    assert power_mod(omega_N2, N1, ell) == 1
    assert power_mod(omega_N1, N2, ell) == 1
    print("Root orders verified ✓")

    # Now verify the full factorization
    errors = 0
    for k1 in range(N1):
        for k2 in range(N2):
            for n1 in range(N1):
                for n2 in range(N2):
                    row = crt_reconstruct(k1, k2)
                    col = ruritanian(n1, n2)
                    actual = V[row][col]
                    expected = (power_mod(omega_N2, k1 * n1, ell) *
                               power_mod(omega_N1, k2 * n2, ell)) % ell
                    if actual != expected:
                        errors += 1

    print(f"\nGood-Thomas factorization errors: {errors} / {D*D}")

    if errors == 0:
        print("\n" + "="*60)
        print("*** GOOD-THOMAS PFA VERIFIED FOR D=96 ***")
        print("="*60)
        print(f"\nV[crt(k1,k2), rur(n1,n2)] = V_32[k1,n1] * V_3[k2,n2]")
        print(f"where V_32[k1,n1] = (ω^3)^{{k1*n1}} mod {ell}")
        print(f"  and V_3[k2,n2] = (ω^32)^{{k2*n2}} mod {ell}")
        print()
        print("=== Implementation as 2-stage MatMul ===")
        print()
        print("Stage 1: 3-point DFT (along k2/n2 component)")
        print(f"  Root: ω^{N1} = {omega_N1} (order 3)")
        print(f"  Matrix V_3 is 3×3, applied to groups of 3 elements")
        print(f"  In slot terms: rotation by 32 positions (σ^32)")
        print(f"  Non-zero diagonals at offsets: 0, 32, 64")
        print(f"  Automorphisms needed: σ^32, σ^64 → 2 auts")
        print()
        print("Stage 2: 32-point DFT (along k1/n1 component)")
        print(f"  Root: ω^{N2} = {omega_N2} (order 32)")
        print(f"  Matrix V_32 is 32×32, applied to groups of 32 elements")
        print(f"  In slot terms: rotation by 3 positions (σ^3)")
        print(f"  Non-zero diagonals at offsets: 0, 3, 6, ..., 93")
        print(f"  With BSGS (baby=6×σ^3, giant=6×σ^18): 11 auts")
        print()
        print("Total automorphisms: 2 + 11 = 13")
        print("Current BSGS for D=96: 18 automorphisms")
        print(f"Reduction: {18-13} fewer = {(18-13)/18*100:.1f}%")
        print()
        print("=== Noise Analysis ===")
        print("Rader97 (failed): 13 sequential stages → 13 key-switch sequences")
        print("Good-Thomas PFA:   2 sequential stages → 2 key-switch sequences")
        print()
        print("Each stage uses standard MatMul1DExec with BSGS/hoisting.")
        print("Stage 1 is extremely sparse (3 diagonals) → minimal noise.")
        print("Stage 2 has 32 diagonals → similar to a standalone 32-point transform.")
        print()
        print("=== Required Permutations ===")
        print("Input permutation (Ruritanian): applied to columns BEFORE transform")
        print("Output permutation (CRT): applied to rows AFTER transform")
        print()
        print("In HElib context:")
        print("  - Permutations are slot rearrangements → can be absorbed into")
        print("    adjacent EvalMap stages or done as a single automorphism step")
        print("  - OR: modify the diagonal constants to account for permutation")
        print("    (fold permutation into the plaintext multipliers)")
        print()

        # Compute the actual diagonal structure for each stage
        print("=== Stage 1 Diagonal Constants ===")
        print("V_3 matrix:")
        for k2 in range(N2):
            row = [power_mod(omega_N1, k2 * n2, ell) for n2 in range(N2)]
            print(f"  row {k2}: {row}")

        print(f"\nStage 2: V_32 is a {N1}×{N1} DFT with root ω^3={omega_N2}")
        print(f"  This is a standard Vandermonde that HElib's MatMul1DExec handles natively.")

    else:
        print(f"\nFACTORIZATION FAILED: {errors} errors")
        # Debug
        for k1 in range(min(2, N1)):
            for k2 in range(N2):
                for n1 in range(min(2, N1)):
                    for n2 in range(N2):
                        row = crt_reconstruct(k1, k2)
                        col = ruritanian(n1, n2)
                        actual = V[row][col]
                        expected = (power_mod(omega_N2, k1*n1, ell) *
                                   power_mod(omega_N1, k2*n2, ell)) % ell
                        if actual != expected:
                            print(f"  k1={k1},k2={k2},n1={n1},n2={n2}: "
                                  f"row={row},col={col}, actual={actual}, expected={expected}")


if __name__ == "__main__":
    main()
