#!/usr/bin/env python3
"""
Key insight from CKKS literature (eprint 2025/1403, Chen-Han 2018, Bossuat et al.):
In CKKS, the CoeffToSlot/SlotToCoeff matrices have MANY zero diagonals because
the DFT matrix for power-of-two has butterfly structure.

For BGV non-power-of-two (m=50731=97×523), the Step2Matrix is a DENSE Vandermonde.
BUT: let's check if the INVERSE Step2Matrix (used for coeffToSlot, the bottleneck)
has any exploitable structure.

Key idea to explore:
1. The Vandermonde matrix V[i,j] = ω^{ij} has the property that V and V^{-1}
   differ only by a scalar: V^{-1} = (1/D) * V^* (conjugate transpose for DFT).
   For exact arithmetic mod ell: V^{-1}[i,j] = (1/D) * ω^{-ij}

2. This means V^{-1} is ALSO a Vandermonde (with ω^{-1} instead of ω).
   So V^{-1} has the SAME diagonal structure as V.

3. The REAL question: can we exploit the CONJUGATE SYMMETRY of the DFT?
   For real-valued inputs: DFT has conjugate symmetry, so half the outputs
   are redundant. This is the "rfft" optimization.

4. For BGV: the inputs to coeffToSlot are NOT real-valued in general.
   BUT: in thin bootstrapping, only SOME slots are used (the "thin" part).
   Can we skip computing unused slots?

5. MORE IMPORTANTLY: Look at the DIAGONAL CONSTANTS of the Step2Matrix.
   In HElib's BSGS, each diagonal d has a constant vector c_d[i] for i=0..D-1.
   If c_d[i] has REPETITIVE structure (e.g., periodic), we can compress the
   computation.

Let me check: for V[i,j] = ω^{ij}, diagonal d has:
  c_d[i] = V[i, (i+d) mod D] = ω^{i * ((i+d) mod D)}

For the INVERSE: V^{-1}[i,j] = (1/D) * ω^{-ij}, diagonal d has:
  c_d[i] = (1/D) * ω^{-i * ((i+d) mod D)}

These are quadratic functions of i in the exponent. No obvious periodicity.

BUT WAIT - there's another approach from the CKKS literature:
"Level-collapsing" (Ma et al. 2024/164): merge consecutive sparse layers
into a single denser layer, trading sparsity for fewer key-switches.

For our case: instead of 2 stages with 2 key-switches, what if we could
find a SINGLE matrix that has fewer than 96 non-zero diagonals?

The Step2Matrix IS dense (96 non-zero diagonals). But what about the
PRODUCT of Step2Matrix with the ADJACENT Step1Matrix or other transforms?

In HElib's thin bootstrapping pipeline:
  slotToCoeff = Step1 (dim=1, D=29) then Step2 (dim=0, D=96)
  coeffToSlot = Step2^{-1} (dim=0, D=96) then Step1^{-1} (dim=1, D=29)

What if we FUSE Step2^{-1} with Step1^{-1}? The fused matrix operates on
BOTH dimensions simultaneously. If the fused matrix has fewer non-zero
"multi-dimensional diagonals", we save key-switches.

Actually, HElib already does this optimization implicitly through the
hypercube structure. Each dimension is handled independently.

Let me think about a completely different angle:

THE REAL INSIGHT FROM 2025/1403:
"By combining the two methods in CoeffsToSlots in a non-trivial manner,
we not only further accelerate the homomorphic linear transformations
and save one level of moduli, but also reduce the total size of evaluation keys."

The key technique in recent CKKS papers is:
1. Split the linear transform into "radix" layers
2. Each layer has few non-zero diagonals (sparse)
3. Use "lazy rescaling" or "level-conserving rescaling" to avoid
   consuming a level between layers

For BGV, we don't have "levels" in the same sense, but we DO have
the analogous concept: each key-switch consumes noise budget.

NEW IDEA: "Lazy key-switching" for the Good-Thomas decomposition.

Instead of doing a full key-switch after Stage 1 (which is what
SparseMatMul1DExec::mul() does via cleanUp()), we can:
1. Apply Stage 1's automorphisms WITHOUT key-switching (just rotate)
2. Accumulate the Stage 1 result as a multi-part ciphertext
3. Apply Stage 2's automorphisms to the multi-part ciphertext
4. Only key-switch ONCE at the very end

This is essentially "hoisting across stages" - the key insight from
the CKKS optimization literature applied to BGV.

In HElib terms: instead of cleanUp() between stages, keep the
ciphertext in "expanded" form (multiple parts) and only reLinearize
at the end.

Let me check if this is feasible by looking at the noise growth.
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

# Check: does the Step2Matrix have any exploitable symmetry?
ell = 97
D = 96
g_ell = find_generator(ell)
omega = g_ell

# Build DFT matrix V[i,j] = omega^{ij} mod ell
V = np.zeros((D, D), dtype=np.int64)
for i in range(D):
    for j in range(D):
        V[i][j] = power_mod(omega, (i*j) % D, ell)

# Check diagonal symmetry: is c_d[i] related to c_{D-d}[i]?
# c_d[i] = V[i, (i+d)%D] = omega^{i*((i+d)%D)}
# c_{D-d}[i] = V[i, (i+D-d)%D] = V[i, (i-d)%D] = omega^{i*((i-d)%D)}

# For the INVERSE matrix: V^{-1}[i,j] = (1/D) * omega^{-ij}
# Diagonal d of V^{-1}: c_d[i] = (1/D) * omega^{-i*((i+d)%D)}

# Key observation: omega^{-1} = omega^{D-1} = omega^95 mod ell
omega_inv = power_mod(omega, D-1, ell)
D_inv = power_mod(D, ell-2, ell)  # D^{-1} mod ell

print(f"omega = {omega}, omega^{{-1}} = {omega_inv}, D^{{-1}} mod {ell} = {D_inv}")

# Check: how many diagonals of V^{-1} are "conjugate pairs"?
# i.e., c_d and c_{D-d} are related by complex conjugation (negating exponent)
# In F_ell, "conjugation" is just omega -> omega^{-1}

# For V^{-1}: diagonal d has c_d[i] = D_inv * omega^{-i*((i+d)%D)}
# diagonal D-d has c_{D-d}[i] = D_inv * omega^{-i*((i-d)%D)}
# = D_inv * omega^{-i^2+id} = D_inv * omega^{id} * omega^{-i^2}
# While c_d[i] = D_inv * omega^{-i^2-id} = D_inv * omega^{-id} * omega^{-i^2}
# So c_{D-d}[i] = c_d[i] * omega^{2id}

# This means: diagonal D-d can be obtained from diagonal d by multiplying
# each entry by omega^{2id}. This is a PLAINTEXT multiplication!
# In HElib terms: if we compute the rotated ciphertext for offset d,
# we can get the contribution of offset D-d by just multiplying by
# a different plaintext constant - NO ADDITIONAL ROTATION NEEDED!

print("\n=== KEY DISCOVERY ===")
print("For the DFT Vandermonde V[i,j] = omega^{ij}:")
print("  Diagonal d: c_d[i] = omega^{i*(i+d) mod D}")
print("  Diagonal D-d: c_{D-d}[i] = omega^{i*(i-d) mod D} = c_d[i] * omega^{-2id}")
print()
print("This means: rotation by d and rotation by D-d share the SAME")
print("automorphism (since sigma^d and sigma^{D-d} = sigma^{-d} are inverses).")
print()
print("Wait - sigma^{D-d} is NOT the same as sigma^d. They are different rotations.")
print("But sigma^{-d} = sigma^{D-d} in Z/DZ.")
print()
print("Actually, the key insight is different:")
print("For a Vandermonde V[i,j] = omega^{ij}, the matrix is SYMMETRIC: V = V^T")
print("(since omega^{ij} = omega^{ji})")
print()

# Verify symmetry
is_symmetric = all(V[i][j] == V[j][i] for i in range(D) for j in range(D))
print(f"V is symmetric: {is_symmetric}")

# For a symmetric matrix, diagonal d and diagonal D-d have a relationship:
# c_d[i] = V[i, (i+d)%D] and c_{D-d}[i] = V[i, (i+D-d)%D] = V[i, (i-d)%D]
# By symmetry V[i,j] = V[j,i]:
# c_{D-d}[i] = V[i, (i-d)%D] = V[(i-d)%D, i] = c_d[(i-d)%D]
# So c_{D-d} is just a CYCLIC SHIFT of c_d!

print("\nFor symmetric Vandermonde:")
print("  c_{D-d}[i] = c_d[(i-d) mod D]")
print("  i.e., diagonal D-d is a cyclic shift of diagonal d!")
print()
print("In BSGS terms: if we have the rotated ciphertext sigma^d(x),")
print("we can compute the contribution of BOTH diagonal d AND diagonal D-d")
print("using the SAME rotated ciphertext, just with different plaintext constants!")
print()
print("This means: we only need rotations for d = 0, 1, 2, ..., D/2 = 48")
print("(not all 96), because each rotation covers TWO diagonals!")
print()
print("=== AUTOMORPHISM REDUCTION ===")
print(f"Standard BSGS for D=96: baby=10, giant=10 → 18 automorphisms")
print(f"With symmetry: only need 48 diagonals → baby=7, giant=7 → 12 automorphisms")
print(f"Reduction: 33%!")
print()

# Verify this concretely
print("Verification: checking c_{D-d}[i] == c_d[(i-d) mod D]...")
errors = 0
for d in range(1, D):
    for i in range(D):
        cd_i = V[i][(i+d) % D]
        cDd_shifted = V[(i-d)%D][((i-d)%D + (D-d)) % D]
        # c_{D-d}[i] = V[i][(i+D-d)%D]
        cDd_i = V[i][(i + D - d) % D]
        # c_d[(i-d)%D] = V[(i-d)%D][((i-d)%D + d) % D] = V[(i-d)%D][i]
        cd_shifted = V[(i-d) % D][i]
        if cDd_i != cd_shifted:
            errors += 1
            if errors <= 3:
                print(f"  MISMATCH at d={d}, i={i}: c_{{D-d}}[{i}]={cDd_i}, c_d[{(i-d)%D}]={cd_shifted}")

if errors == 0:
    print("VERIFIED: c_{D-d}[i] = c_d[(i-d) mod D] for all d, i")
    print()
    print("=== IMPLEMENTATION PLAN ===")
    print("In MatMul1DExec, instead of iterating over all D=96 diagonals,")
    print("iterate over only D/2=48 diagonals (d=0 to 47).")
    print("For each rotation sigma^d(x):")
    print("  1. Multiply by c_d to get contribution of diagonal d")
    print("  2. Multiply by shifted(c_d) to get contribution of diagonal D-d")
    print("  3. Add both to accumulator")
    print()
    print("This halves the number of rotations needed!")
    print("BSGS with 48 effective diagonals: baby=7, giant=7 → 12 auts")
    print("vs current 18 auts → 33% reduction")
    print()
    print("CRITICAL: This is a SINGLE-STAGE optimization (no extra hoisting)!")
    print("It works within the existing MatMul1DExec framework.")
    print("The only change is how diagonal constants are paired.")
else:
    print(f"FAILED: {errors} mismatches")
    print("The simple symmetry doesn't hold. Let me check what DOES hold...")

    # Check the actual relationship
    # For V[i,j] = omega^{ij}: V IS symmetric since ij = ji
    # So V[i,(i+d)%D] = omega^{i*((i+d)%D)}
    # And V[(i-d)%D, i] = omega^{((i-d)%D)*i}
    # These are equal iff i*((i+d)%D) ≡ ((i-d)%D)*i (mod D)
    # i.e., i*(i+d) ≡ (i-d)*i (mod D)
    # i.e., i^2 + id ≡ i^2 - id (mod D)
    # i.e., 2id ≡ 0 (mod D)
    # This is NOT always true! Only when D | 2id.

    print("\nThe relationship c_{D-d}[i] = c_d[(i-d)%D] requires 2id ≡ 0 (mod D)")
    print("This is NOT always true for general i,d.")
    print()
    print("BUT: there's a WEAKER relationship that IS always true:")
    print("sigma^{-d}(x) = sigma^{D-d}(x)")
    print("So rotation by d and rotation by D-d use the SAME key-switch key!")
    print("(sigma^d and sigma^{-d} are handled by the same Galois key in HElib)")
    print()
    print("In HElib, for 'native' dimensions, sigma^d and sigma^{D-d} are")
    print("handled by the SAME rotation. So we already only need D/2 rotations!")
    print()
    print("Let me check if HElib already exploits this...")


if __name__ == "__main__":
    main() if 'main' in dir() else None
