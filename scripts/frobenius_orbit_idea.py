#!/usr/bin/env python3
"""
NEW IDEA: Linearized Polynomial Decomposition for EvalMap

Key mathematical insight from finite field theory:

A "linearized polynomial" (or q-polynomial) over F_q is a polynomial of the form:
  L(X) = a_0 X + a_1 X^q + a_2 X^{q^2} + ... + a_{n-1} X^{q^{n-1}}

These polynomials are F_q-LINEAR maps: L(ax+by) = aL(x) + bL(y) for a,b in F_q.

In our HE setting:
- The slot algebra is F_p[X]/G(X) where p=65537, deg(G)=18
- The Frobenius map σ: x → x^p is an automorphism of the slot algebra
- Applying σ costs exactly 1 key-switch (it's a Galois automorphism)
- The TRACE map Tr(x) = x + x^p + x^{p^2} + ... + x^{p^{17}} is a linearized polynomial
  that maps F_{p^18} → F_p. It costs 17 key-switches naively, but can be computed
  with log2(18) ≈ 5 key-switches using repeated squaring of the Frobenius.

THE KEY OBSERVATION:
The Step2Matrix in HElib computes a linear map over the slot algebra.
Each "diagonal" of the matrix corresponds to a rotation σ^d (along the hypercube dimension).
But there's ANOTHER set of automorphisms available: the Frobenius powers σ_F^k: x → x^{p^k}.

In HElib's framework:
- Rotations along dim=0 (order 96): these are the "expensive" automorphisms used by BSGS
- Frobenius powers (order 18): these are ALSO available as Galois automorphisms

The Step2Matrix entries are elements of F_p[X]/G(X) (degree < 18 polynomials).
When we multiply a ciphertext by a diagonal constant c_d[i], we're doing:
  output_slot[i] += c_d[i] * rotated_slot[i]

where c_d[i] is a polynomial of degree < 18.

NEW IDEA: If c_d[i] can be expressed as a LINEARIZED polynomial in the slot value,
then the multiplication c_d[i] * x can be decomposed as:
  c_d[i] * x = a_0 * x + a_1 * x^p + a_2 * x^{p^2} + ...

where x^{p^k} = Frobenius^k(x) is available via key-switching.

But wait - c_d[i] * x is just a plaintext-ciphertext multiplication, which is ALREADY
cheap (no key-switch needed). So this doesn't help directly.

HOWEVER, there's a deeper insight:

The FULL linear transform L(x) = Σ_d c_d * σ^d(x) can be rewritten as:
  L(x) = Σ_d c_d * σ^d(x)

If we group the diagonals by their Frobenius orbit structure, we might be able
to express L as a composition of:
1. A few rotations (key-switches along dim=0)
2. Frobenius applications (key-switches along the "Frobenius dimension")
3. Plaintext multiplications (free)

The question is: can we find a decomposition that uses FEWER total key-switches
than the standard BSGS?

Let me think about this differently. In HElib's hypercube, there are actually
TWO types of Galois automorphisms:
1. Rotations along dim=0: σ_0^d for d=0,...,95 (order 96)
2. Rotations along dim=1: σ_1^d for d=0,...,28 (order 29)
3. Frobenius: σ_F^k for k=0,...,17 (order 18)

The Step2Matrix for dim=0 uses ONLY type-1 automorphisms.
But what if we could REPLACE some type-1 automorphisms with type-3 (Frobenius)?

The Frobenius σ_F acts on the slot algebra F_{p^18} as x → x^p.
In the "powerful basis" representation, this corresponds to a specific permutation
of the coefficients.

KEY QUESTION: Is there a relationship between the Step2Matrix diagonals and
the Frobenius action that allows us to trade rotations for Frobenius applications?

For the Step2Matrix V[i,j] = (X^{cofactor * g^j})^i mod G(X):
- The entry V[i,j] is a polynomial in F_p[X]/G(X)
- Applying Frobenius to V[i,j]: σ_F(V[i,j]) = V[i,j]^p = (X^{cofactor*g^j*p})^i mod G(X)
  = V[i, j'] where g^{j'} = g^j * p mod ell

So Frobenius permutes the COLUMNS of the Step2Matrix!
Specifically: σ_F maps column j to column (j + ord_p_in_Z96) mod 96.

Since ord_97(p) = 96 and p = 65537, we need: what is the discrete log of p mod 97?
i.e., find k such that g^k ≡ p mod 97.

Let me compute this.
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

def discrete_log(base, target, mod, order):
    """Find k such that base^k ≡ target (mod mod), with base having given order."""
    for k in range(order):
        if power_mod(base, k, mod) == target % mod:
            return k
    return None

ell = 97
p = 65537
D = 96
g = find_generator(ell)
print(f"Generator of (Z/{ell}Z)*: g = {g}")

# Find discrete log of p mod ell w.r.t. generator g
p_mod_ell = p % ell
k_p = discrete_log(g, p_mod_ell, ell, D)
print(f"p mod {ell} = {p_mod_ell}")
print(f"g^{k_p} = {power_mod(g, k_p, ell)} (should be {p_mod_ell})")
print(f"\nFrobenius σ_F permutes column j → column (j + {k_p}) mod {D}")
print(f"Frobenius orbit length on columns: {D // np.gcd(k_p, D)}")

# This means: if we apply Frobenius to the INPUT of the linear transform,
# it's equivalent to cyclically shifting the column index by k_p.
#
# In the diagonal decomposition:
#   L(x) = Σ_d c_d * σ_0^d(x)
#
# Applying Frobenius to x:
#   L(σ_F(x)) = Σ_d c_d * σ_0^d(σ_F(x))
#             = Σ_d c_d * σ_F(σ_0^d(x))  [if σ_0 and σ_F commute]
#             = σ_F(Σ_d σ_F^{-1}(c_d) * σ_0^d(x))
#
# Wait, do σ_0 and σ_F commute? In general, Galois automorphisms of the
# cyclotomic ring commute (the Galois group is abelian for cyclotomic fields).
# So YES, they commute!
#
# This means: σ_F ∘ L = L' ∘ σ_F where L' has diagonals σ_F(c_d).
# Or equivalently: L ∘ σ_F = σ_F ∘ L'' where L'' has diagonals σ_F^{-1}(c_d).
#
# But this doesn't directly help reduce the number of rotations...
#
# UNLESS: some of the diagonal constants c_d are related by Frobenius!
# If c_{d+k_p}[i] = σ_F(c_d[i]) for all i, then we could compute:
#   c_d * σ_0^d(x) + c_{d+k_p} * σ_0^{d+k_p}(x)
#   = c_d * σ_0^d(x) + σ_F(c_d) * σ_0^{d+k_p}(x)
#   = c_d * σ_0^d(x) + σ_F(c_d * σ_0^d(σ_F^{-1}(σ_0^{k_p}(x))))
#
# This is getting complicated. Let me check numerically whether the
# diagonal constants have Frobenius-orbit structure.

print(f"\n=== Checking Frobenius orbit structure of diagonal constants ===")

# Build Step2Matrix over F_ell (scalar approximation)
omega = g
cof = 523 % ell

# Diagonal d has: c_d[i] = V[i, (i+d)%D] = (cof * g^{(i+d)%D})^i mod ell
#                         = cof^i * g^{i*((i+d)%D)} mod ell

# Under Frobenius (which maps column j to column j+k_p):
# The "Frobenius-shifted" diagonal d' = d + k_p should have:
# c_{d+k_p}[i] = cof^i * g^{i*((i+d+k_p)%D)} mod ell
#              = cof^i * g^{i*(i+d)} * g^{i*k_p} mod ell
#              = c_d[i] * g^{i*k_p} mod ell
#              = c_d[i] * (g^{k_p})^i mod ell
#              = c_d[i] * p^i mod ell  (since g^{k_p} = p mod ell)

# So: c_{d+k_p}[i] = c_d[i] * p^i mod ell
# This is a MULTIPLICATIVE relationship between diagonals spaced k_p apart!

print(f"\nTHEOREM: c_{{d+{k_p}}}[i] = c_d[i] * p^i mod {ell}")
print(f"This means: diagonal d+{k_p} = diagonal d * (plaintext vector p^i)")
print(f"\nIn HE terms: if we have σ_0^d(x), we can get the contribution of")
print(f"diagonal d+{k_p} by multiplying by a DIFFERENT plaintext constant,")
print(f"WITHOUT an additional rotation!")
print()

# How many distinct Frobenius orbits are there?
orbit_size = D // np.gcd(k_p, D)
num_orbits = np.gcd(k_p, D)
print(f"Frobenius orbit size: {orbit_size}")
print(f"Number of orbits: {num_orbits}")
print(f"\nThis means: we only need {num_orbits} DISTINCT rotations!")
print(f"Each rotation covers {orbit_size} diagonals (via Frobenius scaling).")
print()

if num_orbits < 18:  # Less than current BSGS
    print(f"*** BREAKTHROUGH: {num_orbits} rotations < 18 (current BSGS)! ***")
    print(f"But we need {orbit_size-1} additional Frobenius applications per orbit.")
    print(f"Total key-switches: {num_orbits} rotations + {num_orbits*(orbit_size-1)} Frobenius")
    print(f"                  = {num_orbits + num_orbits*(orbit_size-1)} total")
    print(f"vs current BSGS: 18 rotations")
else:
    print(f"Number of orbits ({num_orbits}) >= 18, so this doesn't help directly.")
    print(f"But we can combine with BSGS on the orbits.")
    print()
    # With BSGS on num_orbits diagonals:
    import math
    baby = int(math.ceil(math.sqrt(num_orbits)))
    giant = int(math.ceil(num_orbits / baby))
    bsgs_auts = baby + giant - 1
    print(f"BSGS on {num_orbits} orbit representatives: baby={baby}, giant={giant} → {bsgs_auts} auts")
    print(f"Plus Frobenius applications to fill each orbit: {orbit_size-1} per orbit")
    print(f"But Frobenius applications can be done with a SINGLE hoisting!")
    print(f"(All Frobenius powers applied to the same rotated ciphertext)")
    print()
    print(f"Total cost model:")
    print(f"  1 hoisting for rotations: 0.6s")
    print(f"  {bsgs_auts} rotation automorphisms: {bsgs_auts * 0.1:.1f}s")
    print(f"  For each of {num_orbits} rotated ciphertexts:")
    print(f"    {orbit_size} MulAdd operations (one per orbit member): {num_orbits * orbit_size * 0.01:.1f}s")
    print(f"  Total: {0.6 + bsgs_auts*0.1 + num_orbits*orbit_size*0.01:.1f}s")
    print(f"  vs current: {0.6 + 18*0.1 + 96*0.01:.1f}s")

# Verify the relationship numerically
print(f"\n=== Numerical verification ===")
errors = 0
for d in range(D - k_p):
    for i in range(D):
        cd_i = power_mod(cof, i, ell) * power_mod(g, (i * ((i+d) % D)) % D, ell) % ell
        cd_kp_i = power_mod(cof, i, ell) * power_mod(g, (i * ((i+d+k_p) % D)) % D, ell) % ell
        expected = cd_i * power_mod(p_mod_ell, i, ell) % ell
        if cd_kp_i != expected:
            errors += 1
print(f"Verification errors: {errors} (should be 0)")
