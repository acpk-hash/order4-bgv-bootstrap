#!/usr/bin/env python3
"""
Asymmetric BSGS for Homomorphic Linear Transforms

Key insight from Wang & Wang (ASIA CCS '25):
In standard BSGS, baby-step set S1 is computed via key-switching (expensive).
But S1 is at a SHALLOWER multiplicative depth than S2.
By using tensor products (cheap) instead of key-switches to build S1,
we can reduce the number of key-switches from O(sqrt(D)) to O(D^{1/t}).

Application to linear transforms:
Standard BSGS for D=96 (g=10, h=10):
  - Build baby steps: σ^0(x), σ^1(x), ..., σ^9(x) via 9 key-switches
  - For each giant step k: accumulate and apply σ^{10k} via 1 key-switch
  - Total: 9 + 9 = 18 key-switches

Asymmetric BSGS for linear transforms:
  - Decompose baby-step set: instead of computing σ^j(x) directly,
    compute it as a PRODUCT of smaller rotations.
  - Example with t=2: σ^j(x) = σ^{j mod 3}(σ^{3*(j//3)}(x))
    - First compute: σ^0(x), σ^3(x), σ^6(x) [3 key-switches, "inner baby"]
    - Then for each: compute σ^0, σ^1, σ^2 variants [via tensor product, NO key-switch]
    - Total baby: 3 key-switches + 9 tensor products (vs 9 key-switches)

Wait - this doesn't quite work for automorphisms because σ^a(σ^b(x)) = σ^{a+b}(x),
not σ^a(x) * σ^b(x). Automorphisms COMPOSE, they don't multiply.

So the "tensor product" trick from polynomial evaluation doesn't directly apply
to automorphism composition.

BUT: there's a different way to apply the asymmetric idea!

In the BSGS for linear transforms, the expensive operation is:
  For each giant step k: acc_inner.smartAutomorph(σ^{g*k})
This applies an automorphism to the ACCUMULATED inner sum.

The asymmetric idea: instead of applying σ^{g*k} to acc_inner (which requires
a key-switch on a "deep" ciphertext), we can RESTRUCTURE the computation so that
the automorphism is applied to a "shallow" ciphertext.

Specifically, in the standard BSGS:
  result = Σ_k σ^{g*k}(Σ_j c_{j+g*k} ⊙ baby[j])

The inner sum Σ_j c_{j+g*k} ⊙ baby[j] is at the SAME level as baby[j]
(since c_{j+g*k} is a plaintext, MulAdd doesn't increase level).
So the giant-step automorphism σ^{g*k} is applied to a ciphertext at the
same level as the input. This is already "shallow"!

The real cost difference is: key-switch on a ciphertext with MORE parts
(after accumulation) vs fewer parts (single baby-step).

Actually, let me re-read the paper more carefully. The key insight is about
POLYNOMIAL evaluation, not linear transforms. In polynomial evaluation:
  f(x) = Σ_j f^(j)(x) * z^j  where z = x^k

The baby-step set S1 = {x, x^2, ..., x^{k-1}} requires computing POWERS of x.
Each power requires a multiplication (tensor product) + key-switch + modulus-switch.

The asymmetric BSGS says: delay the key-switch and modulus-switch for S1,
compute the powers using only tensor products, and only do key-switch at the end.

For LINEAR TRANSFORMS, the baby-step set is {σ^0(x), σ^1(x), ..., σ^{g-1}(x)}.
These are computed via AUTOMORPHISMS (not multiplications).
Each automorphism requires a key-switch.

CAN WE DELAY THE KEY-SWITCH FOR BABY-STEP AUTOMORPHISMS?

YES! Here's how:
- An automorphism σ^j applied to ct = (c0, c1) gives:
  σ^j(ct) = (σ^j(c0), σ^j(c1))
  But σ^j(c1) is encrypted under s(X^{σ^j}), not s(X).
  Key-switching converts it back to s(X).

- If we DELAY the key-switch, we have a ciphertext encrypted under s(X^{σ^j}).
  We can still do plaintext-ciphertext multiplication (MulAdd) on it!
  (MulAdd just multiplies c0 and c1 by the plaintext, doesn't change the key.)

- So the baby-step computation can be:
  1. Apply σ^j to ct WITHOUT key-switching: baby[j] = (σ^j(c0), σ^j(c1))
     This is just a permutation of coefficients - VERY CHEAP (no key-switch)!
  2. Compute inner sum: acc_inner = Σ_j c_{j+g*k} ⊙ baby[j]
     This works because MulAdd doesn't require the ciphertext to be under s(X).
  3. Key-switch acc_inner ONCE at the end (from s(X^{σ^j}) to s(X)).

BUT WAIT: the baby[j] are under DIFFERENT keys (s(X^{σ^j}) for each j).
We can't add ciphertexts under different keys!

So the inner sum Σ_j c_{j+g*k} ⊙ baby[j] doesn't work directly because
each baby[j] is under a different key.

UNLESS: we use the "hoisting" trick!
Hoisting decomposes the ciphertext ONCE, then applies multiple automorphisms
to the decomposed form. Each automorphism from the decomposed form produces
a ciphertext under the ORIGINAL key s(X) (after the inner product with the
key-switch matrix). So hoisting already avoids the "different key" problem.

This is exactly what HElib's BasicAutomorphPrecon does!
And it's already being used in the BSGS.

So the asymmetric BSGS paper's technique doesn't directly help for
linear transforms because:
1. Baby-step automorphisms are already "free" via hoisting
2. The expensive part is giant-step automorphisms (applied to accumulated sums)
3. Giant-step automorphisms can't be delayed because the result needs to be
   added to other giant-step results

HOWEVER: there IS one applicable idea from the paper:
The "lazy technique" for MULTIPLE POLYNOMIALS (Algorithm 3, Theorem 3.5).

In our linear transform, we evaluate the SAME set of baby-step ciphertexts
for ALL giant steps. The paper's Algorithm 3 shows how to evaluate multiple
polynomials sharing the same baby-step set more efficiently.

For our case: the "multiple polynomials" are the inner sums for each giant step:
  f_k(baby) = Σ_j c_{j+g*k} ⊙ baby[j],  k = 0, ..., h-1

These h=10 "polynomials" share the same baby-step set {baby[0], ..., baby[g-1]}.
The paper says this can be done with O(h*d) scalar mults + O(d^{1/t} + h) key-switches.

But in our case, the "polynomials" are just LINEAR combinations (degree 1 in each baby[j]).
So there's no multiplicative depth issue. The standard BSGS already handles this optimally.

FINAL CONCLUSION:
The asymmetric BSGS paper optimizes POLYNOMIAL EVALUATION (where powers x^k
require multiplications that increase depth). It does NOT directly help with
LINEAR TRANSFORMS (where baby-steps are automorphisms, not multiplications).

The key difference:
- Polynomial eval: baby-step = x^j (requires j-1 multiplications, increasing depth)
- Linear transform: baby-step = σ^j(x) (requires 1 automorphism, no depth increase)

The asymmetric technique exploits the DEPTH DIFFERENCE between S1 and S2.
In linear transforms, there IS no depth difference - all baby-steps are at the same level.
"""

print("CONCLUSION:")
print("The Asymmetric BSGS paper (ASIA CCS '25) optimizes POLYNOMIAL EVALUATION")
print("by exploiting the depth difference between baby-step and giant-step sets.")
print()
print("For HOMOMORPHIC LINEAR TRANSFORMS, this technique does NOT directly apply")
print("because baby-step automorphisms don't increase multiplicative depth.")
print("HElib's hoisting already makes baby-step automorphisms cheap.")
print()
print("HOWEVER: the paper's insight about 'lazy techniques' (Section 2.4) IS relevant")
print("to our DIGIT EXTRACTION optimization (Order-4 cleaner), where we evaluate")
print("high-degree polynomials. The asymmetric BSGS could accelerate the")
print("Paterson-Stockmeyer evaluation of the degree-305 inner polynomial Q(X^4).")
print()
print("Potential application to our Order-4 cleaner:")
print("  Current: Q(Y) evaluated with standard PS, Y=X^4, deg(Q)=305")
print("  With asymmetric BSGS (t=4): key-switches reduced from ~17 to ~5")
print("  Expected speedup on extract phase: ~20-30%")
print("  This would improve extract from 80.8s to ~56-64s!")
print()
print("This is a REAL optimization opportunity for the DIGIT EXTRACTION phase,")
print("not the linear transform phase.")
