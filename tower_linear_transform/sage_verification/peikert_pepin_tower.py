#!/usr/bin/env python3
"""
Peikert-Pepin Tower CRT Transform Implementation (Prototype)

Implements Theorem 5.14 from eprint 2025/1760 for the specific case:
- Number field: K = Q(zeta_97) (97th cyclotomic field)
- Galois group: G = Gal(K/Q) = (Z/97Z)^* ≅ Z/96Z
- Tower: Z/96Z ⊃ <16> (order 6) ⊃ <48> (order 2) ⊃ {0}
  giving degrees [16, 3, 2] at each level
- Plaintext modulus: p = 65537

The CRT transform maps from the "powerful basis" to the "CRT basis"
(evaluation at the 96 roots). The tower decomposition factors this into
stages, each requiring only a small number of automorphisms.

Key formula (Theorem 5.14):
  T_{M/K} = T'_{M/L} * (I_{M/L} ⊗ T_{L/K})

where T_{L/K} is applied to each of the [M:L] blocks independently,
and T'_{M/L} is a "lifted" transform with diagonal blocks of size [M:L].
"""

import numpy as np
from typing import List, Tuple

def mod_inv(a: int, p: int) -> int:
    """Modular inverse using Fermat's little theorem."""
    return pow(int(a) % int(p), int(p) - 2, int(p))

def factorize(n: int) -> List[int]:
    """Return prime factors of n."""
    factors = []
    d = 2
    while d * d <= n:
        while n % d == 0:
            factors.append(d)
            n //= d
        d += 1
    if n > 1:
        factors.append(n)
    return factors

def primitive_root_mod(n: int) -> int:
    """Find smallest primitive root modulo prime n."""
    fs = set(factorize(n - 1))
    for g in range(2, n):
        if all(pow(g, (n - 1) // int(f), n) != 1 for f in fs):
            return g
    raise ValueError(f"No primitive root found for {n}")

def find_working_field(ell: int, D: int) -> int:
    """Find a prime q where both ell-th and D-th roots of unity exist."""
    from math import gcd
    lcm_val = ell * D // gcd(ell, D)
    for k in range(2, 1000000):
        q = lcm_val * k + 1
        if q < 2**31 and all(q % d != 0 for d in range(2, min(int(q**0.5) + 1, 100000))):
            return q
    raise ValueError("No suitable working field found")


class CyclotomicCRTTransform:
    """
    Implements the CRT transform for Q(zeta_ell) with tower decomposition.

    For ell=97: G = (Z/97Z)^* is cyclic of order 96.
    We choose a subgroup tower based on the factorization 96 = d1 * d2 * ...
    """

    def __init__(self, ell: int = 97, p: int = 65537, cofactor: int = 523):
        self.ell = ell
        self.p = p
        self.D = ell - 1  # = 96
        self.cofactor = cofactor
        self.cofactor_mod = cofactor % ell

        # Find working field
        self.q = find_working_field(ell, self.D)
        print(f"Working field: F_{self.q}")

        # Primitive roots
        self.gen_q = primitive_root_mod(self.q)
        self.zeta = pow(self.gen_q, (self.q - 1) // ell, self.q)

        # Generator of (Z/97Z)^*
        self.g = primitive_root_mod(ell)  # = 5
        self.g_inv = mod_inv(self.g, ell)

        # HElib's reps: reps[i] = g_inv^i mod ell
        self.reps = [pow(int(self.g_inv), i, ell) for i in range(self.D)]

        # Evaluation points
        self.points = [pow(self.zeta, (self.reps[j] * self.cofactor_mod) % ell, self.q)
                      for j in range(self.D)]

        # Build the full Step2Matrix (Vandermonde)
        self.A = self._build_vandermonde()

    def _build_vandermonde(self) -> np.ndarray:
        """Build the 96x96 Vandermonde matrix."""
        A = np.zeros((self.D, self.D), dtype=np.int64)
        for j in range(self.D):
            power = 1
            for i in range(self.D):
                A[i][j] = power
                power = (power * self.points[j]) % self.q
        return A

    def _mat_mul_mod(self, A: np.ndarray, B: np.ndarray) -> np.ndarray:
        """Matrix multiplication modulo q."""
        n = A.shape[0]
        C = np.zeros((n, n), dtype=np.int64)
        for i in range(n):
            for j in range(n):
                val = 0
                for k in range(n):
                    val = (val + int(A[i][k]) * int(B[k][j])) % self.q
                C[i][j] = val
        return C

    def _mat_inv_mod(self, M: np.ndarray) -> np.ndarray:
        """Compute matrix inverse modulo q using Gauss-Jordan."""
        n = M.shape[0]
        aug = np.zeros((n, 2*n), dtype=np.int64)
        for i in range(n):
            for j in range(n):
                aug[i][j] = int(M[i][j]) % self.q
            aug[i][n+i] = 1

        for col in range(n):
            pivot = -1
            for row in range(col, n):
                if aug[row][col] % self.q != 0:
                    pivot = row
                    break
            if pivot == -1:
                raise ValueError("Matrix is singular")
            aug[[col, pivot]] = aug[[pivot, col]]
            inv_pivot = mod_inv(int(aug[col][col]), self.q)
            aug[col] = (aug[col] * inv_pivot) % self.q
            for row in range(n):
                if row != col and aug[row][col] % self.q != 0:
                    factor = int(aug[row][col])
                    aug[row] = (aug[row] - factor * aug[col]) % self.q

        return aug[:, n:] % self.q

    def _count_diagonals(self, M: np.ndarray) -> Tuple[int, set]:
        """Count non-zero diagonals of a matrix."""
        n = M.shape[0]
        diags = set()
        for i in range(n):
            for j in range(n):
                if int(M[i][j]) % self.q != 0:
                    diags.add((j - i) % n)
        return len(diags), diags

    def compute_tower_factorization(self, tower: List[int]):
        """
        Compute the tower factorization of the Step2Matrix.

        tower = [d1, d2, ...] where D = d1 * d2 * ...

        The factorization is: A = M_last * ... * M_2 * M_1
        where M_k has non-zero diagonals at offsets that are multiples of (d1*...*d_{k-1}).

        Stage k uses offsets: {0, s, 2s, ..., (d_k-1)*s} where s = d1*...*d_{k-1}
        """
        D = self.D
        assert np.prod(tower) == D, f"Tower {tower} doesn't multiply to {D}"

        print(f"\nTower decomposition: {tower}")
        print(f"D = {' × '.join(map(str, tower))} = {D}")

        # Compute the stage offsets
        stages = []
        stride = 1
        for k, dk in enumerate(tower):
            offsets = [stride * j for j in range(dk)]
            stages.append(offsets)
            print(f"  Stage {k+1}: {dk} diagonals at offsets {offsets[:8]}{'...' if dk > 8 else ''}")
            stride *= dk

        # Verify coverage: all D offsets should be reachable
        from itertools import product as cartprod
        all_offsets = set()
        for combo in cartprod(*stages):
            all_offsets.add(sum(combo) % D)
        assert len(all_offsets) == D, f"Only {len(all_offsets)} offsets covered (need {D})"
        print(f"  Coverage: {len(all_offsets)}/{D} offsets ✓")

        # Now compute the actual stage matrices.
        # Strategy: compute stages from RIGHT to LEFT.
        # M_1 is computed first, then M_2 = A * M_1^{-1} restricted to stage-2 offsets, etc.
        #
        # For a 2-stage tower [d1, d2]:
        # A = M2 * M1
        # M1 has diags at {0, 1, ..., d1-1}
        # M2 has diags at {0, d1, 2*d1, ..., (d2-1)*d1}
        #
        # We CHOOSE M1 and then M2 = A * M1^{-1}.
        # The question: what M1 gives M2 with the right sparsity?
        #
        # KEY INSIGHT from Peikert-Pepin:
        # M1 = T_{L/K} ⊗ I_{[M:L]} (the "lower" CRT transform tensored with identity)
        # This means M1 is BLOCK DIAGONAL with [M:L] = d2 identical blocks of size d1.
        # Each block is the d1-point CRT transform of the "lower" subfield.

        # For our case with tower [d1, d2]:
        # M1 is block-diagonal: d2 blocks of size d1, each block = T_{d1-point CRT}
        # The d1-point CRT is a d1×d1 Vandermonde at d1 specific points.

        # These d1 points are the "coarse" evaluation points for the subfield.
        # Specifically: they are the d1 distinct values of zeta^(g^{d2*k}) for k=0,...,d1-1
        # (the elements of the subgroup of order d1 in (Z/ell)^*).

        if len(tower) == 2:
            return self._compute_2stage(tower[0], tower[1])
        elif len(tower) == 3:
            return self._compute_3stage(tower)
        else:
            raise NotImplementedError(f"Tower of length {len(tower)} not implemented")

    def _compute_2stage(self, d1: int, d2: int):
        """Compute 2-stage tower factorization [d1, d2]."""
        D = self.D
        q = self.q

        print(f"\n  Computing 2-stage factorization [{d1}, {d2}]...")

        # M1 is block-diagonal: d2 blocks of size d1.
        # Each block is a d1×d1 Vandermonde at the "coarse" points.
        #
        # The coarse points for block b (b=0,...,d2-1) are:
        # points in positions {b, b+d2, b+2*d2, ..., b+(d1-1)*d2}
        # i.e., the points at indices that are ≡ b (mod d2).
        #
        # Actually, for the Peikert-Pepin construction:
        # M1 = I_{d2} ⊗ V_{d1} where V_{d1} is a d1×d1 sub-Vandermonde.
        # This means: M1[i][j] ≠ 0 only if (j-i) mod D < d1 AND i//d1 == j//d1
        # Wait, that's not right for the diagonal representation.

        # Let me think about this differently.
        # In the "diagonal offset" representation:
        # M1 has non-zero diags at offsets {0, 1, ..., d1-1}
        # M1[i][(i+d) mod D] = c_d[i] for d = 0,...,d1-1
        #
        # For M1 = I_{d2} ⊗ V_{d1} (block diagonal with d2 blocks of V_{d1}):
        # Block b contains rows {b*d1, b*d1+1, ..., b*d1+d1-1}
        # Within block b: M1[b*d1+r][b*d1+c] = V_{d1}[r][c]
        # In offset terms: offset = (b*d1+c) - (b*d1+r) = c - r (mod D)
        # Since 0 ≤ r,c < d1, offset ∈ {-(d1-1), ..., d1-1} mod D = {0,...,d1-1, D-d1+1,...,D-1}
        # That's 2*d1-1 possible offsets, not d1!

        # Hmm, this doesn't match. The issue is that "block diagonal" in the
        # NATURAL ordering doesn't correspond to "few diagonals" in the CYCLIC ordering.

        # The correct interpretation for the CYCLIC (slot rotation) setting:
        # M1 with diags at {0,1,...,d1-1} means:
        # M1[i][(i+d) mod D] = c_d[i] for d=0,...,d1-1
        # This is a "d1-tap FIR filter" in the cyclic slot domain.

        # For this to equal I_{d2} ⊗ V_{d1}, we need the slot ordering to align
        # with the block structure. This requires a SPECIFIC slot ordering.

        # In HElib, slots are ordered by the generator g_s of the dimension.
        # Slot i corresponds to the automorphism sigma^i (rotation by i).
        # The subgroup H1 = <d1> = {0, d1, 2*d1, ...} corresponds to rotations
        # by multiples of d1.

        # For M1 to be "block diagonal" in the rotation sense, it must commute
        # with rotations by d1. A matrix that commutes with rot^{d1} has
        # non-zero diags only at multiples of... no, that's the OPPOSITE.
        # A matrix with diags at {0,...,d1-1} does NOT commute with rot^{d1}.

        # I think the correct construction is:
        # M1 has diags at {0, d2, 2*d2, ..., (d1-1)*d2}  (stride d2, not stride 1!)
        # M2 has diags at {0, 1, 2, ..., d2-1}
        # Then A = M2 * M1 covers all offsets {a*d2 + b : a=0..d1-1, b=0..d2-1} = all D.

        # Let me try BOTH orderings and see which one works.

        # Option A: M1 diags at {0,1,...,d1-1}, M2 diags at {0,d1,...,(d2-1)*d1}
        # Option B: M1 diags at {0,d2,...,(d1-1)*d2}, M2 diags at {0,1,...,d2-1}

        # For Option B with d1=6, d2=16:
        # M1 diags at {0, 16, 32, 48, 64, 80} — stride 16
        # M2 diags at {0, 1, 2, ..., 15}

        # For Option A with d1=6, d2=16:
        # M1 diags at {0, 1, 2, 3, 4, 5}
        # M2 diags at {0, 6, 12, 18, ..., 90}

        # Let me try Option B first (M1 has large-stride diags).
        # M1 with diags at {0, 16, 32, 48, 64, 80} commutes with rot^16
        # (since all offsets are multiples of 16... no, multiples of 16 means
        # it commutes with rot^1 restricted to the quotient Z/96Z / <16>...
        # actually a matrix with diags at multiples of d2=16 DOES commute with rot^1
        # on the "coarse" level. Let me just try numerically.

        # NUMERICAL APPROACH: Try M1 = sub-Vandermonde at coarse points
        # Coarse points: take every d2-th point
        # points[0], points[d2], points[2*d2], ..., points[(d1-1)*d2]

        # Build M1 as block-circulant with d1 blocks, each of size d2
        # M1[i][(i + k*d2) mod D] = coarse_point[i mod d2]^k for k=0,...,d1-1

        # Actually let me just try the simplest thing:
        # Set M1[i][(i+k*d2) mod D] = points[i]^k for k=0,...,d1-1
        # This makes M1 a "strided Vandermonde" in the rotation domain.

        M1 = np.zeros((D, D), dtype=np.int64)
        for i in range(D):
            for k in range(d1):
                j = (i + k * d2) % D
                # M1[i][j] = points[i]^k (the i-th point raised to power k)
                M1[i][j] = pow(int(self.points[i]), k, q)

        # Check if M1 is invertible
        try:
            M1_inv = self._mat_inv_mod(M1)
        except ValueError:
            print("  M1 (strided Vandermonde) is SINGULAR. Trying different construction...")
            return None

        # Compute M2 = A * M1^{-1}
        M2 = self._mat_mul_mod(self.A, M1_inv)

        # Check diagonal structure of M2
        n_diags_M2, diags_M2 = self._count_diagonals(M2)

        # We want M2 to have diags only at {0, 1, ..., d2-1}
        wanted_diags = set(range(d2))
        unwanted = diags_M2 - wanted_diags

        print(f"  M1: diags at multiples of {d2} (stride {d2}), {d1} diagonals")
        print(f"  M2 = A * M1^{{-1}}: {n_diags_M2} non-zero diagonals")
        print(f"  Wanted ({{0,...,{d2-1}}}): {len(wanted_diags & diags_M2)} present")
        print(f"  Unwanted: {len(unwanted)} extra diagonals")

        if len(unwanted) == 0:
            print(f"  ✓ FACTORIZATION WORKS! A = M2 * M1")
            print(f"  M1: {d1} diags (offsets {list(range(0, D, d2))[:6]}...)")
            print(f"  M2: {d2} diags (offsets {list(range(d2))})")

            # Verify: M2 * M1 == A
            product = self._mat_mul_mod(M2, M1)
            match = all((product[i][j] - self.A[i][j]) % q == 0
                       for i in range(D) for j in range(D))
            print(f"  Verification M2*M1 == A: {match}")
            return M1, M2
        else:
            print(f"  ✗ Factorization failed. M2 has extra diagonals.")
            # Try the other direction
            print(f"\n  Trying reversed: M1 diags at {{0,...,{d1-1}}}, M2 diags at multiples of {d1}...")

            M1b = np.zeros((D, D), dtype=np.int64)
            for i in range(D):
                for k in range(d1):
                    j = (i + k) % D
                    M1b[i][j] = pow(int(self.points[i]), k, q)

            try:
                M1b_inv = self._mat_inv_mod(M1b)
            except ValueError:
                print("  M1b is SINGULAR.")
                return None

            M2b = self._mat_mul_mod(self.A, M1b_inv)
            n_diags_M2b, diags_M2b = self._count_diagonals(M2b)
            wanted_b = set(range(0, D, d1))
            unwanted_b = diags_M2b - wanted_b

            print(f"  M2b: {n_diags_M2b} non-zero diagonals")
            print(f"  Wanted (multiples of {d1}): {len(wanted_b & diags_M2b)} present")
            print(f"  Unwanted: {len(unwanted_b)} extra diagonals")

            if len(unwanted_b) == 0:
                print(f"  ✓ REVERSED FACTORIZATION WORKS!")
                return M1b, M2b
            else:
                print(f"  ✗ Both directions failed.")
                return None


def main():
    print("=" * 70)
    print("Peikert-Pepin Tower CRT Transform - Prototype Implementation")
    print("=" * 70)

    transform = CyclotomicCRTTransform(ell=97, p=65537, cofactor=523)

    # Try various tower decompositions
    # Tower [6, 16]: 2 stages
    result = transform.compute_tower_factorization([6, 16])

    if result is None:
        print("\n\nTrying tower [16, 6]...")
        result = transform.compute_tower_factorization([16, 6])

    if result is None:
        print("\n\nTrying tower [3, 32]...")
        result = transform.compute_tower_factorization([3, 32])

    if result is None:
        print("\n\nTrying tower [32, 3]...")
        result = transform.compute_tower_factorization([32, 3])


if __name__ == "__main__":
    main()
