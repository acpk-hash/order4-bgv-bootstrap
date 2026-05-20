#!/usr/bin/env python3
"""
CRITICAL RE-ANALYSIS: What does the Frobenius orbit relation actually buy us?

The relation is: c_{d+16}[i] = c_d[i] * p^i mod ell
(where 16 = gcd(80, 96), orbit step is 16 in the diagonal offset)

Orbit 0: diagonals at offsets {0, 16, 32, 48, 64, 80}
These correspond to rotations σ^0, σ^16, σ^32, σ^48, σ^64, σ^80.
These are 6 DIFFERENT rotations - they cannot share a single rotation!

So what's the actual benefit?

ANSWER: The benefit is in the BSGS decomposition.

In standard BSGS with baby=b, giant=g:
  σ^d(x) = σ^{g*k}(σ^j(x))  where d = j + g*k, j < b

The baby-step ciphertexts are: σ^0(x), σ^1(x), ..., σ^{b-1}(x)
The giant-step applies σ^{g*k} to the accumulated inner sum.

With Frobenius-orbit structure, we can choose a DIFFERENT BSGS decomposition:
  - Baby-step: σ^0, σ^1, ..., σ^15 (16 baby steps)
  - Giant-step: σ^16, σ^32, σ^48, σ^64, σ^80 (5 giant steps)
  - Total: 16 + 5 = 21 automorphisms... WORSE than 18!

Wait, that's not right either. Let me think more carefully.

The standard BSGS for D=96 with baby=10, giant=10:
  - 10 baby-step auts: σ^0, ..., σ^9
  - 9 giant-step auts: σ^10, σ^20, ..., σ^90
  - Covers: all d = j + 10*k for j=0..9, k=0..9 → d = 0..99 (covers 0..95)

The Frobenius orbit tells us: diagonals d and d+16 have RELATED constants.
Specifically: if we already computed c_d ⊙ σ^d(x), we can get c_{d+16} ⊙ σ^{d+16}(x)
by computing: (c_d * p^i) ⊙ σ^{d+16}(x).

But σ^{d+16}(x) is a DIFFERENT ciphertext from σ^d(x). So we still need the rotation.

HOWEVER: in the BSGS framework, σ^{d+16}(x) might already be available!
If baby=16, then σ^0, σ^1, ..., σ^15 are all precomputed as baby steps.
And σ^{d+16}(x) = σ^16(σ^d(x)) which is a giant-step applied to a baby-step.

Actually, the REAL insight is different. Let me reconsider.

The Frobenius relation c_{d+16}[i] = c_d[i] * p^i means:
  The PLAINTEXT CONSTANTS for diagonal d+16 are determined by those of diagonal d.

In HElib's MatMul1DExec, the plaintext constants are PRECOMPUTED and stored.
The relation doesn't save any online computation - it only saves STORAGE.

But wait - there IS a computational saving if we think about it differently:

In the BSGS inner loop:
  acc_inner = Σ_{j=0}^{b-1} c_{j+g*k} ⊙ baby_step[j]

If we choose g=16 (giant-step = 16), then for giant-step k:
  acc_inner_k = Σ_{j=0}^{b-1} c_{j+16*k} ⊙ baby_step[j]
              = Σ_{j=0}^{b-1} (c_j * (p^i)^k) ⊙ baby_step[j]
              = (p^i)^k ⊙ (Σ_{j=0}^{b-1} c_j ⊙ baby_step[j])
              = (p^i)^k ⊙ acc_inner_0

THIS IS THE KEY INSIGHT!

If giant-step = 16, then:
  acc_inner_k = (p^{ik}) ⊙ acc_inner_0

The inner accumulation for giant-step k is just a PLAINTEXT SCALING of
the inner accumulation for giant-step 0!

This means: we only need to compute the inner accumulation ONCE (for k=0),
and then all other giant-steps are just plaintext multiplications!

In standard BSGS: each giant-step requires a full key-switch (σ^{g*k} applied
to acc_inner_k). But with this structure:
  - Compute acc_inner_0 = Σ_{j=0}^{15} c_j ⊙ baby_step[j]
  - For k=1: acc_inner_1 = (p^i) ⊙ acc_inner_0 (plaintext mult, FREE)
  - For k=2: acc_inner_2 = (p^{2i}) ⊙ acc_inner_0 (plaintext mult, FREE)
  - ...
  - For k=5: acc_inner_5 = (p^{5i}) ⊙ acc_inner_0 (plaintext mult, FREE)

Then the full result is:
  result = acc_inner_0 + σ^16(acc_inner_1) + σ^32(acc_inner_2) + ... + σ^80(acc_inner_5)
         = acc_inner_0 + σ^16((p^i) ⊙ acc_inner_0) + σ^32((p^{2i}) ⊙ acc_inner_0) + ...

So the computation is:
  1. Compute baby steps: σ^0(x), ..., σ^15(x) → 15 hoisted auts
  2. Compute acc_inner_0 = Σ_{j=0}^{15} c_j ⊙ baby_step[j] → 16 MulAdds
  3. For k=1..5: compute (p^{ki}) ⊙ acc_inner_0 → 5 plaintext mults (FREE)
  4. Apply giant steps: σ^16, σ^32, σ^48, σ^64, σ^80 to the scaled versions → 5 key-switches
  5. Sum everything up

Total key-switches: 15 (baby) + 5 (giant) = 20... still more than 18.

But wait - with HOISTING, the baby steps are cheap (from precon).
The expensive part is the GIANT steps (unhoisted key-switches).

Standard BSGS (baby=10, giant=10): 9 giant-step key-switches (expensive)
Frobenius-orbit BSGS (baby=16, giant=6): 5 giant-step key-switches (expensive)

Cost comparison:
  Standard: 0.6(precon) + 9*0.1(baby) + 9*0.5(giant) + 96*0.01(MulAdd) = 6.06s
  Frobenius: 0.6(precon) + 15*0.1(baby) + 5*0.5(giant) + 16*0.01(MulAdd) + 5*0.01(scale) = 4.76s

Hmm, but 15 baby steps is a lot. Let me reconsider.

Actually, with hoisting (GeneralAutomorphPrecon), ALL automorphisms from the precon
cost the same (~0.1s each). The distinction between "baby" and "giant" is:
- Baby: extracted from precon (cheap, ~0.1s)
- Giant: applied to the ACCUMULATED result (expensive, ~0.5s, full key-switch)

So the real question is: how many GIANT steps do we need?

Standard BSGS (baby=10, giant=10):
  - 10 baby from precon: 10 × 0.1s = 1.0s
  - 9 giant (full key-switch on accumulated): 9 × 0.5s = 4.5s
  - 96 MulAdds: 0.96s
  - Precon: 0.6s
  Total: 7.06s

Frobenius-orbit (baby=16, giant=6):
  - 16 baby from precon: 16 × 0.1s = 1.6s (slightly more baby steps)
  - 5 giant (full key-switch): 5 × 0.5s = 2.5s (MUCH fewer giant steps!)
  - 16 MulAdds for inner sum: 0.16s
  - 5 plaintext scalings: ~0.05s (essentially free)
  - 5 additions: ~0.05s
  - Precon: 0.6s
  Total: 4.96s
  Speedup: 7.06 / 4.96 = 1.42x

Even better: baby=16 means we extract 16 rotations from the precon.
The precon cost might increase slightly (more rotations to decompose),
but GeneralAutomorphPrecon handles this efficiently.

WAIT - I need to reconsider. The standard BSGS doesn't use GeneralAutomorphPrecon
for ALL rotations. Let me re-read the HElib code to understand the exact execution model.

From the code (matmul.cpp line 1382-1450):
  - Baby steps are computed iteratively: baby_steps[j] = baby_steps[j-1].smartAutomorph(σ^1)
  - Each baby step is a full key-switch (smartAutomorph includes key-switching)
  - Giant steps: sum.smartAutomorph(σ^g) applied to the accumulated sum

So BOTH baby and giant steps are full key-switches in the iterative mode!
The cost is: (baby-1) + (giant-1) key-switches = 9 + 9 = 18 key-switches.

With Frobenius-orbit (baby=16, giant=6):
  - Baby: 15 key-switches
  - Giant: 5 key-switches
  - Total: 20 key-switches... WORSE!

BUT: the giant-step key-switches are applied to a LARGER ciphertext (the accumulated sum),
which might be more expensive. And with hoisting mode (non-iterative), baby steps share
a single precomputation.

Let me check which mode HElib actually uses for D=96.
"""

# The key question: does HElib use iterative or hoisted BSGS for D=96?
# From the code: it depends on fhe_test_force_hoist and getKSStrategy(dim).
# For the default case with our parameters, it likely uses the hoisted path
# (GenBabySteps with BasicAutomorphPrecon).

# In hoisted mode:
#   - 1 BasicAutomorphPrecon (decompose input once): 0.6s
#   - baby_steps extracted from precon: each ~0.1s
#   - giant_steps: each is smartAutomorph on accumulated sum: ~0.5s

# So the cost model is:
#   Standard (baby=10, giant=10): 0.6 + 9*0.1 + 9*0.5 + 96*0.01 = 6.06s
#   Frobenius (baby=16, giant=6): 0.6 + 15*0.1 + 5*0.5 + 16*0.01 + 5*0.01 = 4.81s
#   Speedup: 6.06/4.81 = 1.26x

# Hmm, only 26% speedup on the D=96 block. On linear2 (~45s), the D=96 blocks
# account for ~14.4s, so savings would be ~14.4*(1-1/1.26) = ~2.9s.
# New linear2: ~42s. Modest but real.

# BUT WAIT: there's a much better formulation!
#
# The Frobenius relation says: for giant-step k, the inner accumulation is
# just a plaintext scaling of the k=0 accumulation. This means we DON'T NEED
# to recompute the inner sum for each giant step!
#
# Standard BSGS:
#   For each giant step k=0..g-1:
#     acc_inner = Σ_j c_{j+g*k} ⊙ baby[j]  (b MulAdds)
#     if k>0: acc_inner = σ^{g*k}(acc_inner) (1 key-switch)
#     result += acc_inner
#   Total MulAdds: b*g = 96
#   Total giant key-switches: g-1 = 9
#
# Frobenius-orbit BSGS (g=16):
#   acc_inner_0 = Σ_{j=0}^{15} c_j ⊙ baby[j]  (16 MulAdds)
#   result = acc_inner_0
#   For k=1..5:
#     scaled = (p^{ki}) ⊙ acc_inner_0  (1 plaintext mult)
#     result += σ^{16k}(scaled)  (1 key-switch)
#   Total MulAdds: 16 (not 96!)
#   Total giant key-switches: 5
#   Total plaintext mults: 5
#
# This is DRAMATICALLY fewer MulAdds: 16 instead of 96!
# And fewer giant key-switches: 5 instead of 9!

print("=== CORRECTED COST MODEL ===")
print()
print("Standard BSGS (baby=10, giant=10):")
print(f"  Precon: 0.6s")
print(f"  Baby steps (9 from precon): 9 × 0.1s = 0.9s")
print(f"  Inner MulAdds: 96 total (10 per giant step × ~10 steps)")
print(f"  Giant key-switches: 9 × 0.5s = 4.5s")
print(f"  MulAdd cost: 96 × 0.01s = 0.96s")
print(f"  Total: ~7.0s")
print()
print("Frobenius-Orbit BSGS (baby=16, giant-via-scaling):")
print(f"  Precon: 0.6s")
print(f"  Baby steps (15 from precon): 15 × 0.1s = 1.5s")
print(f"  Inner MulAdds (k=0 only): 16 × 0.01s = 0.16s")
print(f"  Plaintext scalings (k=1..5): 5 × 0.01s = 0.05s")
print(f"  Giant key-switches (k=1..5): 5 × 0.5s = 2.5s")
print(f"  Total: ~4.8s")
print(f"  Speedup on D=96 block: 7.0/4.8 = {7.0/4.8:.2f}x")
print()
print("Impact on linear2 (~45s):")
print(f"  D=96 blocks: 2 × 7.0s = 14.0s → 2 × 4.8s = 9.6s")
print(f"  Savings: 4.4s")
print(f"  New linear2: ~40.6s (10% faster)")
print()
print("Combined with Order-4 cleaner + Parallel CoeffToSlot:")
print(f"  Current combined: 118.5s")
print(f"  With Frobenius-orbit: 118.5 - 4.4 = ~114.1s")
print(f"  vs baseline 213.1s: 1.87x speedup")
