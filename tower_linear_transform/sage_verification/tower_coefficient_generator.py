#!/usr/bin/env python3
"""
Tower Decomposition Coefficient Generator.

Builds Step2Matrix for ell=97, p=65537, then factors it as M2 × M1 where:
- M1: 6 diagonals at offsets {0, 16, 32, 48, 64, 80} (Frobenius subgroup)
- M2: 16 diagonals at offsets {0, 1, ..., 15} (coset representatives)

The Step2Matrix is a 96×96 matrix over F_p where:
  S[i][j] = eval_point[j]^i mod p
and eval_point[j] = g^(cofactor * reps[j]) mod ell, evaluated as a scalar
via the embedding F_p[X]/(factor_j(X)) -> F_p (choosing a root).

For the "thin" Step2Matrix used in HElib, each slot sees a SCALAR matrix.
The 96 positions correspond to the 96 automorphisms σ^0, ..., σ^95.
The matrix entry S[i][j] is the j-th evaluation point raised to power i.
"""

import json
import numpy as np
from pathlib import Path


def mod_inv(a, p):
    return pow(int(a) % int(p), int(p) - 2, int(p))


def build_step2_matrix(ell=97, p=65537):
    """
    Build the 96×96 Step2Matrix over F_p.

    The evaluation points are the 96-th roots of unity in F_p
    (since 96 | p-1 = 65536, all 96-th roots exist in F_p).

    Specifically: point[j] = omega^j where omega is a primitive 96-th root of unity mod p.
    The Step2Matrix is the Vandermonde: S[i][j] = point[j]^i mod p.
    """
    D = ell - 1  # 96

    # Find a primitive 96-th root of unity mod p
    # Since p-1 = 65536 = 2^16, and 96 = 2^5 * 3... wait, 3 does not divide p-1!
    # So 96-th roots of unity do NOT all exist in F_p.
    # Only (gcd(96, p-1) = gcd(96, 65536) = 32)-th roots exist.

    # This means the Step2Matrix is NOT a simple scalar Vandermonde over F_p.
    # It's a Vandermonde over the EXTENSION field F_{p^6} (since ord_97(p) = 6).

    # For HElib's "thin" bootstrapping, the Step2Matrix operates on the
    # POLYNOMIAL representation. Each "slot" contains a degree-5 polynomial
    # (element of F_{p^6}), and the matrix maps between coefficient and
    # evaluation representations.

    # The correct construction: the evaluation points are elements of F_{p^6},
    # represented as polynomials mod an irreducible factor of Phi_97(X).

    # For our purposes, we work in F_p[X]/(Phi_97(X)) and the Step2Matrix
    # entries are POLYNOMIALS (vectors of length 96 over F_p).

    # However, for the DIAGONAL representation used in HElib's MatMul1D,
    # what matters is the diagonal values. The diagonal at offset k is:
    # diag_k = [S[0][(0+k)%96], S[1][(1+k)%96], ..., S[95][(95+k)%96]]

    # For a Vandermonde S[i][j] = point[j]^i:
    # diag_k[i] = point[(i+k)%96]^i

    # The "points" are the images of ζ under the 96 automorphisms:
    # point[j] = σ^j(ζ) = ζ^{g^j} where g=5 is the generator of (Z/97Z)*

    # In F_p[X]/(Phi_97(X)): ζ = X, so point[j] = X^{g^j mod 97} mod Phi_97(X)
    # These are polynomials, not scalars!

    # For the tower decomposition, we need to work with the FULL polynomial
    # representation. This makes the computation more complex.

    # ALTERNATIVE APPROACH: Work numerically in a large enough field.
    # Use F_q where q = p^6 and 97 | q-1 (so 97-th roots of unity exist).
    # Then the Step2Matrix becomes a scalar Vandermonde over F_q.

    # For VERIFICATION purposes, we can work mod a small prime q where
    # 97 | q-1 and check the factorization structure.

    # Let's find such a prime: need q ≡ 1 mod 97 and q prime.
    # q = 97*k + 1 for some k. Try: 97*2+1=195 (not prime), 97*4+1=389 (prime!)
    # But 389 is too small. Let's use q = 97*674+1 = 65379... not prime.
    # Actually, let's just use the polynomial representation directly.

    pass  # Will implement below


def poly_mul_mod_phi97(f, g, p):
    """Multiply polynomials f*g mod (Phi_97(X), p). Degree bound: 95."""
    D = 96
    # Phi_97(X) = X^96 + X^95 + ... + X + 1
    result = [0] * (2 * D - 1)
    for i in range(min(len(f), D)):
        if f[i] == 0:
            continue
        for j in range(min(len(g), D)):
            if g[j] == 0:
                continue
            result[i + j] = (result[i + j] + f[i] * g[j]) % p

    # Reduce mod Phi_97: X^96 ≡ -(1 + X + ... + X^95) mod Phi_97
    for i in range(len(result) - 1, D - 1, -1):
        if result[i] != 0:
            c = result[i]
            for j in range(D):
                result[i - D + j] = (result[i - D + j] - c) % p
            result[i] = 0

    return [x % p for x in result[:D]]


def poly_pow_mod_phi97(f, exp, p):
    """Compute f^exp mod (Phi_97(X), p)."""
    D = 96
    result = [0] * D
    result[0] = 1
    base = list(f[:D]) + [0] * (D - len(f))

    while exp > 0:
        if exp & 1:
            result = poly_mul_mod_phi97(result, base, p)
        base = poly_mul_mod_phi97(base, base, p)
        exp >>= 1

    return result


def build_step2_poly_matrix(p=65537, ell=97):
    """
    Build Step2Matrix where each entry is a polynomial in F_p[X]/(Phi_97).

    S[i][j] = (X^{g^j mod 97})^i mod Phi_97(X) mod p

    where g=5 is the primitive root of (Z/97Z)*.

    Returns: 96×96 matrix, each entry is a list of 96 coefficients.
    """
    D = ell - 1  # 96
    g = 5  # primitive root mod 97

    # Compute evaluation points: point[j] = X^{g^j mod 97} mod Phi_97
    # These are monomials X^{exp} where exp = g^j mod 97
    exponents = []
    gj = 1
    for j in range(D):
        exponents.append(gj)
        gj = (gj * g) % ell

    # point[j] is the polynomial with a single 1 at position exponents[j]
    # But X^97 ≡ 1 mod Phi_97 (since Phi_97 | X^97 - 1), so X^{exp} for exp < 97
    # is just the monomial X^{exp} (no reduction needed since exp < 96...
    # wait, exp can be up to 96, and Phi_97 has degree 96.
    # X^96 ≡ -(1+X+...+X^95) mod Phi_97

    points = []
    for j in range(D):
        exp = exponents[j]
        poly = [0] * D
        if exp < D:  # exp < 96
            poly[exp] = 1
        else:  # exp = 96 (only possible if g^j mod 97 = 96... but 96 < 97 so exp <= 96)
            # X^96 mod Phi_97 = -(1+X+...+X^95)
            poly = [(p - 1)] * D  # all -1 mod p
        points.append(poly)

    print(f"Built {D} evaluation points (polynomials mod Phi_97)")

    # Build Step2Matrix: S[i][j] = point[j]^i mod Phi_97
    # This is expensive: 96*96 polynomial exponentiations
    print("Building Step2Matrix (96×96 polynomial entries)...")

    S = [[None] * D for _ in range(D)]

    for j in range(D):
        # Compute point[j]^i for i=0,...,95
        power = [0] * D
        power[0] = 1  # point[j]^0 = 1
        S[0][j] = list(power)

        for i in range(1, D):
            power = poly_mul_mod_phi97(power, points[j], p)
            S[i][j] = list(power)

        if (j + 1) % 16 == 0:
            print(f"  Column {j+1}/{D} done")

    return S, exponents


def extract_scalar_diagonals(S, D=96, p=65537):
    """
    For the thin Step2Matrix, extract the scalar diagonal representation.

    In HElib's thin bootstrapping, the matrix multiplication is:
    result_slot[i] = Σ_k diag[k][i] * input_slot[(i+k) % D]

    where diag[k][i] = S[i][(i+k) % D] (a polynomial, but for thin case
    we take the constant term or the trace).

    Actually for thin bootstrapping, each slot has degree d=6, so the
    matrix is 16×16 (nslots/d = 96/6 = 16)...

    No wait - the Step2Matrix for dimension 0 (the 97-component) has size D=96.
    It operates on the 96 coefficients of the polynomial representation.
    Each diagonal entry is a SCALAR (element of F_p), not a polynomial.

    The scalar diagonal at offset k:
    diag_k[i] = S[i][(i+k)%D] evaluated as a scalar.

    But S[i][j] is a polynomial! For the scalar case, we need the
    "constant coefficient" or some specific evaluation.

    Actually, I think the issue is that for the 97-component with d=6,
    the Step2Matrix is 16×16 (operating on the 16 slots), not 96×96.
    Let me reconsider.
    """
    # TODO: Need to understand the exact HElib matrix structure better
    pass


def verify_tower_structure(p=65537, ell=97):
    """
    Verify the tower decomposition structure numerically.

    Strategy: Work in a field where 97-th roots of unity are scalars.
    Find prime q such that 97 | q-1, then the Step2Matrix is a scalar
    Vandermonde over F_q.
    """
    D = ell - 1  # 96
    g = 5  # primitive root mod 97

    # Find a working prime q with 97 | q-1
    # We need q large enough and q ≡ 1 mod 97
    q = 97 * 2 + 1  # 195, not prime
    for k in range(2, 10000):
        q = 97 * k + 1
        if is_prime(q):
            break

    print(f"Working prime: q = {q} (q-1 = {q-1}, 97 | q-1: {(q-1)%97 == 0})")

    # Find a primitive 97-th root of unity mod q
    # omega = g_q^{(q-1)/97} where g_q is a primitive root mod q
    g_q = find_primitive_root_fast(q)
    omega = pow(g_q, (q - 1) // 97, q)
    assert pow(omega, 97, q) == 1 and pow(omega, 1, q) != 1
    print(f"Primitive 97-th root of unity: omega = {omega}")

    # Evaluation points: point[j] = omega^{g^j mod 97}
    exponents = []
    gj = 1
    for j in range(D):
        exponents.append(gj)
        gj = (gj * g) % ell

    points = [pow(omega, exp, q) for exp in exponents]

    # Build scalar Vandermonde: S[i][j] = points[j]^i mod q
    S = np.zeros((D, D), dtype=np.int64)
    for j in range(D):
        power = 1
        for i in range(D):
            S[i][j] = power
            power = (power * points[j]) % q

    print(f"Built {D}×{D} Vandermonde over F_{q}")

    # Frobenius subgroup: {0, 16, 32, 48, 64, 80}
    frob_offsets = [0, 16, 32, 48, 64, 80]
    coset_offsets = list(range(16))

    # Extract M1: block-diagonal with 16 blocks of 6×6
    # M1 has diagonals at offsets {0, 16, 32, 48, 64, 80}
    # M1[i][j] ≠ 0 only if (j-i)%96 ∈ frob_offsets, i.e., j%16 == i%16

    # Strategy: M1 is the "intra-coset" Vandermonde.
    # Within coset c (positions {c, c+16, c+32, c+48, c+64, c+80}),
    # M1 restricted to this coset is a 6×6 Vandermonde with points
    # {points[c], points[c+16], points[c+32], points[c+48], points[c+64], points[c+80]}

    # Build M1 as block-diagonal
    M1 = np.zeros((D, D), dtype=np.int64)
    for c in range(16):  # for each coset
        coset_indices = [c + 16 * k for k in range(6)]
        coset_points = [points[idx] for idx in coset_indices]
        # 6×6 Vandermonde within this coset
        for local_i in range(6):
            for local_j in range(6):
                global_i = coset_indices[local_i]
                global_j = coset_indices[local_j]
                M1[global_i][global_j] = pow(coset_points[local_j], local_i, q)

    # Verify M1 has support only on offsets {0, 16, 32, 48, 64, 80}
    m1_offsets = set()
    for i in range(D):
        for j in range(D):
            if M1[i][j] != 0:
                m1_offsets.add((j - i) % D)
    print(f"M1 diagonal offsets: {sorted(m1_offsets)}")
    assert m1_offsets == set(frob_offsets), f"M1 offsets mismatch! Got {sorted(m1_offsets)}"
    print("  ✓ M1 has correct diagonal structure")

    # Compute M2 = S × M1^{-1} mod q
    # First invert M1 mod q
    M1_inv = matrix_inv_mod(M1, q)
    if M1_inv is None:
        print("  ✗ M1 is not invertible mod q!")
        return False

    # M2 = S × M1^{-1}
    M2 = matrix_mul_mod(S, M1_inv, q)

    # Check M2 diagonal offsets
    m2_offsets = set()
    for i in range(D):
        for j in range(D):
            if M2[i][j] % q != 0:
                m2_offsets.add((j - i) % D)
    print(f"M2 diagonal offsets: {sorted(m2_offsets)}")

    if m2_offsets.issubset(set(coset_offsets)):
        print("  ✓ M2 has correct diagonal structure (offsets ⊆ {0,...,15})")
    else:
        extra = m2_offsets - set(coset_offsets)
        print(f"  ✗ M2 has extra offsets: {sorted(extra)}")
        return False

    # Verify: M2 × M1 == S
    product = matrix_mul_mod(M2, M1, q)
    match = all(product[i][j] % q == S[i][j] % q for i in range(D) for j in range(D))
    print(f"Verification M2 × M1 == S: {'✓' if match else '✗'}")

    if match:
        # Extract diagonal values
        m1_diags = {}
        for offset in frob_offsets:
            diag = [int(M1[i][(i + offset) % D]) for i in range(D)]
            m1_diags[offset] = diag

        m2_diags = {}
        for offset in coset_offsets:
            diag = [int(M2[i][(i + offset) % D]) % q for i in range(D)]
            m2_diags[offset] = diag

        print(f"\nExtracted M1 diagonals: {len(m1_diags)} (offsets {sorted(m1_diags.keys())})")
        print(f"Extracted M2 diagonals: {len(m2_diags)} (offsets {sorted(m2_diags.keys())})")

        # Save to JSON
        output = {
            "q": q,
            "ell": ell,
            "D": D,
            "frob_offsets": frob_offsets,
            "coset_offsets": coset_offsets,
            "m1_diags": {str(k): v for k, v in m1_diags.items()},
            "m2_diags": {str(k): v for k, v in m2_diags.items()},
        }

        out_path = Path(__file__).parent.parent / "tower_coeffs_test.json"
        with open(out_path, "w") as f:
            json.dump(output, f, indent=2)
        print(f"\nSaved test coefficients to {out_path}")

        return True

    return False


def is_prime(n):
    if n < 2:
        return False
    if n < 4:
        return True
    if n % 2 == 0 or n % 3 == 0:
        return False
    i = 5
    while i * i <= n:
        if n % i == 0 or n % (i + 2) == 0:
            return False
        i += 6
    return True


def find_primitive_root_fast(p):
    """Find smallest primitive root mod prime p."""
    phi = p - 1
    factors = set()
    n = phi
    d = 2
    while d * d <= n:
        while n % d == 0:
            factors.add(d)
            n //= d
        d += 1
    if n > 1:
        factors.add(n)

    for g in range(2, p):
        if all(pow(g, phi // f, p) != 1 for f in factors):
            return g
    return None


def matrix_inv_mod(M, q):
    """Invert matrix M mod q using Gaussian elimination."""
    D = len(M)
    # Augment [M | I]
    aug = [[0] * (2 * D) for _ in range(D)]
    for i in range(D):
        for j in range(D):
            aug[i][j] = int(M[i][j]) % q
        aug[i][D + i] = 1

    # Forward elimination
    for col in range(D):
        # Find pivot
        pivot = -1
        for row in range(col, D):
            if aug[row][col] % q != 0:
                pivot = row
                break
        if pivot == -1:
            return None  # singular

        # Swap
        aug[col], aug[pivot] = aug[pivot], aug[col]

        # Scale pivot row
        inv_pivot = pow(aug[col][col] % q, q - 2, q)
        for j in range(2 * D):
            aug[col][j] = (aug[col][j] * inv_pivot) % q

        # Eliminate
        for row in range(D):
            if row == col:
                continue
            factor = aug[row][col] % q
            if factor == 0:
                continue
            for j in range(2 * D):
                aug[row][j] = (aug[row][j] - factor * aug[col][j]) % q

    # Extract inverse
    result = np.zeros((D, D), dtype=np.int64)
    for i in range(D):
        for j in range(D):
            result[i][j] = aug[i][D + j] % q

    return result


def matrix_mul_mod(A, B, q):
    """Multiply matrices A × B mod q."""
    D = len(A)
    C = np.zeros((D, D), dtype=np.int64)
    for i in range(D):
        for j in range(D):
            s = 0
            for k in range(D):
                s += int(A[i][k]) * int(B[k][j])
            C[i][j] = s % q
    return C


if __name__ == "__main__":
    print("=" * 60)
    print("TOWER DECOMPOSITION COEFFICIENT GENERATOR")
    print("=" * 60)
    print()
    verify_tower_structure()
