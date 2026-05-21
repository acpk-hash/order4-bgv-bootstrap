#!/usr/bin/env python3
"""
New idea: Fused Single-Stage Good-Thomas BSGS

Instead of 2 sequential stages (which doubles the hoisting cost),
implement a SINGLE MatMul1DExec that exploits the Good-Thomas structure
WITHIN its BSGS execution.

Key insight: In standard BSGS for D=96 (baby=10, giant=10):
  - 96 non-zero diagonals → 96 plaintext-ciphertext multiplications
  - 10 baby-step automorphisms + 9 giant-step automorphisms = 19 total

In Good-Thomas-aware BSGS:
  - The 96 diagonals have ALGEBRAIC RELATIONSHIPS between them
  - Specifically: diagonal d and diagonal d+32 differ only by a factor of ω_3
    (because of the 3-point DFT structure)
  - This means: instead of 96 independent MulAdd operations, we can do
    32 groups of 3, where each group shares computation

But wait - this doesn't reduce automorphisms, only plaintext multiplications.
And plaintext-ciphertext mults are already cheap (~7ms each).

Let me think about what ACTUALLY dominates the cost:
- From the timer: mul_MatMul1DExec = 7.36s per call
- keySwitchPart = 0.628s per call, 83 total = 52.1s
- automorph (in matmul) = 0.105s per call, 151 total = 15.9s
- BasicAutomorphPrecon = 0.601s per call, 8 total = 4.8s

The 83 keySwitchPart calls correspond to:
- 6 MatMul1DExec calls (for the 2 dimensions × 3 repeats... no, 6 calls for 1 repeat)
- Each call does: baby-step hoisting (1 precon) + giant-step automorphisms

For D=96 with baby=10, giant=10:
- 1 hoisting precomputation (BasicAutomorphPrecon) = 0.6s
- 10 baby-step automorphisms (from precon) = 10 × 0.1s = 1.0s
- 9 giant-step smartAutomorphs = 9 × 0.5s = 4.5s (these include key-switching)
- 96 MulAdd operations = 96 × 0.01s = ~1s
Total per D=96 call: ~7.1s ✓ (matches observed 7.36s)

The giant-step automorphisms are the expensive part (4.5s out of 7.36s).
Each giant-step does a full key-switch (0.5s).

With Good-Thomas (baby=3 for 3-point, then baby=6 for 32-point):
If we could do it in ONE stage with shared hoisting:
- 1 hoisting precomputation = 0.6s
- Baby-steps for ALL needed rotations: σ^1, σ^2, σ^3, ..., σ^{12} = 12 from precon = 1.2s
- Giant-steps: σ^{13}, σ^{26}, ..., = fewer needed

Actually, the real question is: what's the OPTIMAL baby/giant split for D=96
when we know the matrix is sparse (only 3+32=35 non-zero diagonals in the
factored form)?

Wait - the FULL Step2Matrix has 96 non-zero diagonals (it's dense).
The factored form has 3 + 32 = 35 diagonals across 2 stages.
But in a single stage, it's still 96 diagonals.

The ONLY way to reduce automorphisms in a single stage is to have fewer
non-zero diagonals. And the Step2Matrix IS dense (96 non-zero diagonals).

So the fundamental tradeoff is:
- 1 stage, 96 diagonals: 18 auts, 1 hoisting = 0.6 + 1.0 + 4.5 + 1.0 = 7.1s
- 2 stages, 3+32 diagonals: 2+11=13 auts, 2 hoistings = 1.2 + 1.3 + 3.0 + 0.4 = 5.9s

Hmm, the 2-stage should actually be faster! Let me check why it wasn't.

The issue might be that SparseMatMul1DExec doesn't use hoisting (it uses
buildGeneralAutomorphPrecon which is different from BasicAutomorphPrecon).

Let me check the actual execution path for SparseMatMul1DExec.
"""

# The key numbers:
# BasicAutomorphPrecon (hoisting): 0.6s per call
# automorph from precon: 0.1s per call
# smartAutomorph (giant-step, includes key-switch): 0.5s per call
# MulAdd: ~0.01s per call

# Standard BSGS for D=96 (baby=10, giant=10):
# 1 hoisting + 10 baby auts + 9 giant auts + 96 MulAdds
# = 0.6 + 10*0.1 + 9*0.5 + 96*0.01 = 0.6 + 1.0 + 4.5 + 0.96 = 7.06s

# Good-Thomas 2-stage:
# Stage 1 (3 diagonals): 1 hoisting + 3 baby auts + 0 giant + 3*96 MulAdds...
# Wait, SparseMatMul1DExec uses GeneralAutomorphPrecon, not BasicAutomorphPrecon.
# GeneralAutomorphPrecon does hoisting for ALL offsets at once.
# For 3 offsets: 1 hoisting + 3 automorphs from precon + 3*96/3 MulAdds
# Actually: for D=96 with 3 diagonals at {0, 32, 64}:
#   1 precon + 2 automorphs (offset 32, 64; offset 0 is identity) + 96 MulAdds
#   = 0.6 + 2*0.1 + 96*0.01 = 0.6 + 0.2 + 0.96 = 1.76s

# Stage 2 (32 diagonals at {0,3,6,...,93}):
#   1 precon + 31 automorphs + 96*32/96 MulAdds...
#   Wait, SparseMatMul1DExec iterates over offsets, not baby/giant.
#   For each offset, it calls precon->automorph(offset).
#   With GeneralAutomorphPrecon, this is hoisting-based: 1 decomposition, then
#   each automorph is cheap (just a rotation of the decomposed form).
#   So: 1 precon(0.6s) + 31 automorphs(31*0.1s) + 32*96 MulAdds... no.
#   Actually each MulAdd is per-position, so it's 32 MulAdds total (one per diagonal).
#   Wait no - MulAdd multiplies a constant vector by a ciphertext and adds to accumulator.
#   There are 32 non-zero diagonals, so 32 MulAdd calls.
#   = 0.6 + 31*0.1 + 32*0.01 = 0.6 + 3.1 + 0.32 = 4.02s

# Total 2-stage: 1.76 + 4.02 = 5.78s
# vs baseline: 7.06s
# Expected speedup: 18%

# But observed: Good-Thomas took 46.36s for linear2 vs baseline 45.27s (SLOWER!)
# This means my cost model is wrong, or there's additional overhead.

# Let me check: the Good-Thomas run has 8 BasicAutomorphPrecon calls (same as baseline!)
# and 151 automorphs in matmul (same as baseline 453/3=151 per repeat).
# This suggests Good-Thomas is NOT being triggered for the D=96 dimension!

# Wait - looking more carefully:
# Baseline (repeat=3): mul_MatMul1DExec: 128.705 / 18 = 7.15s
# Good-Thomas (repeat=1): mul_MatMul1DExec: 44.1725 / 6 = 7.36s
#
# 6 calls means: 2 dimensions × 3 stages? No, for thin boot there are
# nfactors=2 dimensions (D=96 and D=29). With Good-Thomas on D=96,
# the D=96 dimension should use ThinStep2GoodThomasExec (not MatMul1DExec).
# But we see 6 MatMul1DExec calls...
#
# Possible explanation: Good-Thomas is only enabled for forward (not invert),
# and the 6 calls are: 3 forward + 3 inverse? No, repeat=1.
# Actually for thin boot: slotToCoeff (1 call) + coeffToSlot (1 call) = 2 calls
# Each has 2 dimensions = 4 MatMul1DExec calls.
# Plus the Good-Thomas stages use SparseMatMul1DExec (not MatMul1DExec).
#
# So the 6 MatMul1DExec calls are the NON-Good-Thomas dimensions.
# The Good-Thomas stages would show up as SparseMatMul1DExec.
# But there's no SparseMatMul1DExec timer in the output!
#
# This means: Good-Thomas is NOT being triggered. The condition check
# might be failing (e.g., the 'inflate' parameter).

print("DIAGNOSIS: Good-Thomas may not be triggering for the D=96 dimension.")
print("Check: the 'inflate=true' case in the constructor.")
print("The buildThinStep2GoodThomasExec returns nullptr when inflate=true!")
print()
print("Looking at the code:")
print("  if (!thinStep2GoodThomasEnabled() || invert || inflate)")
print("    return nullptr;")
print()
print("The FIRST call (sz==nfactors case) passes inflate=true!")
print("The SECOND loop passes inflate=false.")
print()
print("For m=50731=97×523, nfactors=2:")
print("  dim=1 (D=96): called with inflate=true in the first block")
print("  dim=0 (D=29): called with inflate=false in the loop")
print()
print("So Good-Thomas is only attempted for D=29 (which fails the coprime check)")
print("and NOT for D=96 (because inflate=true)!")
print()
print("FIX: Remove the 'inflate' restriction from buildThinStep2GoodThomasExec,")
print("or handle the inflate case inside the constructor.")
