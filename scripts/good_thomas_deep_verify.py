#!/usr/bin/env python3
"""
Deeper analysis of Step2Matrix structure for D=96.

The Step2Matrix is V[i,j] = point_j^i mod p, where:
  point_j = (cofactor * g^j) mod ell, for j=0,...,95
  g is a generator of (Z/97Z)*, so g has order 96 mod 97.

Key insight: the COLUMN index j corresponds to the group element g^j.
Under CRT Z/96Z ≅ Z/32Z × Z/3Z, we write j = CRT(a, b).
Then g^j = g^{CRT(a,b)} = (g^{D2_inv_mod_D1 * D2})^a * (g^{D1_inv_mod_D2 * D1})^b

This means point_{CRT(a,b)} = point_0^{...} which has multiplicative structure.

The correct factorization is NOT on the matrix indices directly,
but on the EXPONENTS in the Vandermonde.

Let's think about it differently:
  V[i, j] = ω^{i*j} where ω = point_1 = g mod ell (primitive 96th root in F_ell)

Wait - that's only true if point_j = ω^j. Let me check.
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
    """Find a primitive root mod prime m."""
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
        is_gen = True
        for f in factors:
            if power_mod(g, phi // f, m) == 1:
                is_gen = False
                break
        if is_gen:
            return g
    return None


def main():
    p = 65537
    ell = 97
    D = 96
    D1 = 32
    D2 = 3
    cofactor = 523  # m / ell = 50731 / 97

    g_ell = find_generator(ell)
    print(f"Generator of (Z/{ell}Z)*: g = {g_ell}")
    print(f"Order of g: {D}")

    # reps[j] = g^j mod ell
    reps = [power_mod(g_ell, j, ell) for j in range(D)]

    # points[j] = (cofactor * reps[j]) mod ell
    cof_mod_ell = cofactor % ell
    print(f"cofactor mod ell = {cof_mod_ell}")
    points = [(cof_mod_ell * reps[j]) % ell for j in range(D)]

    # Key observation: points[j] = (cof * g^j) mod ell
    # So points[j] = points[0] * (g^j / g^0) = ... no, it's multiplicative:
    # points[j] = cof * g^j mod ell
    # points[j+1] = cof * g^{j+1} = points[j] * g mod ell

    # The Vandermonde matrix is:
    # A[i,j] = points[j]^i = (cof * g^j)^i = cof^i * g^{ij} mod ell
    #
    # So A[i,j] = cof^i * g^{ij} mod ell
    #
    # The factor cof^i only depends on the ROW, not the column!
    # Therefore: A = diag(1, cof, cof^2, ..., cof^{D-1}) * V
    # where V[i,j] = g^{ij} mod ell is a PURE DFT matrix!

    print("\n=== KEY INSIGHT ===")
    print("A[i,j] = (cof * g^j)^i = cof^i * g^{ij} mod ell")
    print("A = diag(cof^0, cof^1, ..., cof^{D-1}) * DFT_matrix")
    print("The diagonal scaling does NOT affect the automorphism structure!")
    print("(It's just a plaintext-ciphertext multiplication per row)")
    print()

    # Verify this
    # Build pure DFT: V[i,j] = g^{ij} mod ell
    V_dft = np.zeros((D, D), dtype=np.int64)
    for i in range(D):
        for j in range(D):
            V_dft[i][j] = power_mod(g_ell, i * j, ell)

    # Build diagonal: D_mat[i] = cof^i mod ell
    diag_cof = [power_mod(cof_mod_ell, i, ell) for i in range(D)]

    # Verify A = diag(cof^i) * V_dft
    errors = 0
    for i in range(D):
        for j in range(D):
            expected = (diag_cof[i] * V_dft[i][j]) % ell
            actual = (power_mod(cof_mod_ell * power_mod(g_ell, j, ell), i, ell))
            if expected != actual:
                errors += 1
    print(f"Verification A = diag(cof^i) * DFT: errors = {errors}")

    # Now: DFT matrix V[i,j] = g^{ij} mod ell, where g has order 96 mod 97.
    # Under CRT: i = 32*c + 3*d (mod 96), j = 32*a + 3*b (mod 96)
    # (using CRT isomorphism Z/96Z ≅ Z/32Z × Z/3Z)
    #
    # V[i,j] = g^{ij} = g^{(32c+3d)(32a+3b)} = g^{1024ca + 96da + 96cb + 9db}
    #         = g^{1024ca} * g^{96(da+cb)} * g^{9db}
    #         = g^{1024ca} * 1 * g^{9db}    (since g^96 = 1)
    #         = (g^{1024})^{ca} * (g^9)^{db}
    #
    # g^{1024 mod 96} = g^{1024 - 10*96} = g^{64}  (order 96/gcd(64,96) = 96/32 = 3)
    # g^9 has order 96/gcd(9,96) = 96/3 = 32
    #
    # So V[CRT(c,d), CRT(a,b)] = (g^64)^{ca} * (g^9)^{db}
    #                            = V_3[c,a] * V_32[d,b]
    # where V_3[c,a] = (g^64)^{ca} is a 3×3 DFT
    # and   V_32[d,b] = (g^9)^{db} is a 32×32 DFT

    print("\n=== Good-Thomas Factorization of DFT ===")
    g64 = power_mod(g_ell, 64, ell)
    g9 = power_mod(g_ell, 9, ell)
    print(f"g^64 mod {ell} = {g64}, order = {96 // np.gcd(64, 96)} (should be 3)")
    print(f"g^9 mod {ell} = {g9}, order = {96 // np.gcd(9, 96)} (should be 32)")

    # Verify order
    assert power_mod(g64, 3, ell) == 1, f"g^64 doesn't have order 3!"
    assert power_mod(g64, 1, ell) != 1, f"g^64 has order 1!"
    assert power_mod(g9, 32, ell) == 1, f"g^9 doesn't have order 32!"
    assert power_mod(g9, 16, ell) != 1, f"g^9 doesn't have order 16!"
    print("Orders verified!")

    # Build CRT mapping
    # CRT: (c, d) -> c * D2 * inv(D2, D1) + d * D1 * inv(D1, D2) mod D
    inv_D2_D1 = pow(D2, -1, D1)  # inv(3, 32) = 11
    inv_D1_D2 = pow(D1, -1, D2)  # inv(32, 3) = inv(2, 3) = 2
    print(f"\nCRT parameters: inv({D2},{D1})={inv_D2_D1}, inv({D1},{D2})={inv_D1_D2}")

    def to_crt(c, d):
        return (c * D2 * inv_D2_D1 + d * D1 * inv_D1_D2) % D

    def from_crt(idx):
        return idx % D1, idx % D2  # (c mod 32, d mod 3)

    # Verify factorization: V[to_crt(c,d), to_crt(a,b)] = (g^64)^{ca} * (g^9)^{db} mod ell
    print("\nVerifying V[CRT(c,d), CRT(a,b)] = V_3[c,a] * V_32[d,b]...")
    errors = 0
    for c in range(D2):  # 0,1,2
        for d in range(D1):  # 0,...,31
            for a in range(D2):
                for b in range(D1):
                    i = to_crt(c, d)  # Note: CRT maps (c,d) -> row index
                    j = to_crt(a, b)  # CRT maps (a,b) -> col index
                    # Wait - need to be careful about which component maps to which
                    # i mod D1 = c*D2*inv_D2_D1 mod D1... this is getting confused.
                    # Let me use the direct formula:
                    # V[i,j] = g^{ij}, and i*j mod 96 determines the value.
                    actual = V_dft[i][j]
                    # Under CRT: i = to_crt(c,d), j = to_crt(a,b)
                    # i*j mod 96: we need to compute this
                    # The factorization says: g^{ij} = g^{64*ca} * g^{9*db} when
                    # i ≡ c (mod 3), i ≡ d (mod 32), j ≡ a (mod 3), j ≡ b (mod 32)
                    # Actually the correct mapping is:
                    # i mod 3 = "3-component of i", i mod 32 = "32-component of i"
                    i_mod3 = i % D2
                    i_mod32 = i % D1
                    j_mod3 = j % D2
                    j_mod32 = j % D1
                    expected = (power_mod(g64, i_mod3 * j_mod3, ell) *
                               power_mod(g9, i_mod32 * j_mod32, ell)) % ell
                    if actual != expected:
                        errors += 1

    print(f"Factorization errors: {errors} / {D*D}")

    if errors == 0:
        print("\n*** PERFECT FACTORIZATION VERIFIED! ***")
        print(f"\nV_DFT[i,j] = g^{{ij}} = (g^64)^{{(i mod 3)(j mod 3)}} * (g^9)^{{(i mod 32)(j mod 32)}} mod {ell}")
        print(f"\nThis means the 96-point DFT factors as:")
        print(f"  DFT_96 = (DFT_3 ⊗ I_32) · (I_3 ⊗ DFT_32)")
        print(f"  under the natural CRT decomposition of indices.")
        print(f"\nFor the full Step2Matrix:")
        print(f"  Step2Matrix = diag(cof^i) · DFT_96")
        print(f"              = diag(cof^i) · (DFT_3 ⊗ I_32) · (I_3 ⊗ DFT_32)")
        print(f"\nThe diagonal can be absorbed into Stage 1 or Stage 2.")
        print(f"\n=== Automorphism Count ===")
        print(f"Stage 1 (DFT_3 along mod-3 component):")
        print(f"  Automorphisms needed: σ^32, σ^64 (rotations by 32 and 64)")
        print(f"  Count: 2")
        print(f"\nStage 2 (DFT_32 along mod-32 component):")
        print(f"  Automorphisms needed: σ^3, σ^6, ..., σ^93")
        print(f"  With BSGS (baby=6 of σ^3, giant=6 of σ^18): 6+6-1 = 11")
        print(f"  Count: 11")
        print(f"\nTotal: 2 + 11 = 13 automorphisms")
        print(f"Current BSGS: 18 automorphisms")
        print(f"Reduction: 28%")
        print(f"\nCRITICAL: Only 2 sequential MatMul stages (vs Rader97's 13)!")
        print(f"Each stage inherits HElib's BSGS + hoisting.")
        print(f"Noise: ~2 key-switch sequences (vs 13 for Rader97).")
    else:
        print(f"\nFactorization failed with {errors} errors. Investigating...")
        # Print a few mismatches
        count = 0
        for i in range(D):
            for j in range(D):
                i_mod3 = i % D2
                i_mod32 = i % D1
                j_mod3 = j % D2
                j_mod32 = j % D1
                expected = (power_mod(g64, i_mod3 * j_mod3, ell) *
                           power_mod(g9, i_mod32 * j_mod32, ell)) % ell
                actual = V_dft[i][j]
                if actual != expected and count < 5:
                    print(f"  i={i}, j={j}: actual g^{(i*j)%96}={actual}, "
                          f"expected g64^{i_mod3*j_mod3}*g9^{i_mod32*j_mod32}={expected}")
                    print(f"    i*j mod 96 = {(i*j)%96}")
                    print(f"    (i%3)*(j%3) = {i_mod3*j_mod3}, (i%32)*(j%32) = {i_mod32*j_mod32}")
                    print(f"    64*(i%3)*(j%3) + 9*(i%32)*(j%32) mod 96 = {(64*i_mod3*j_mod3 + 9*i_mod32*j_mod32)%96}")
                    print(f"    Should equal i*j mod 96 = {(i*j)%96}")
                    count += 1


if __name__ == "__main__":
    main()
