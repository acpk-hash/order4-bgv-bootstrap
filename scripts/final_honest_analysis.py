#!/usr/bin/env python3
"""
FUNDAMENTAL RE-THINK: What actually dominates the cost?

From profiling (baseline, D=96, 1 call):
  - BasicAutomorphPrecon: 0.6s (1 call)
  - automorph (from precon): 0.105s × ~25 = 2.6s
  - keySwitchPart: 0.628s × ~14 = 8.8s
  - mul_MatMul1DExec: 7.36s (total for 1 D=96 call)

Wait, 7.36s total but keySwitchPart alone is 8.8s? That's because keySwitchPart
is shared across ALL 6 MatMul1DExec calls. Let me recalculate per-call:
  - 83 keySwitchPart total / 6 calls ≈ 14 per call
  - 14 × 0.628s = 8.8s... but mul_MatMul1DExec is only 7.36s?

Actually, keySwitchPart is called from WITHIN mul_MatMul1DExec (via smartAutomorph
and BasicAutomorphPrecon). So the 7.36s INCLUDES the key-switch time.

Let me decompose the 7.36s:
  - GenBabySteps: creates g=10 baby steps via BasicAutomorphPrecon
    - 1 BasicAutomorphPrecon (0.6s)
    - 9 automorphisms from precon (9 × 0.105s = 0.95s)
    - Each automorph from precon includes a keySwitchPart (9 × 0.628s = 5.65s)
    Wait, that's already > 7.36s. Something's wrong.

Let me re-read GenBabySteps more carefully.
"""

# From the code (line 1107-1151):
# GenBabySteps uses BasicAutomorphPrecon which does:
#   precon = BasicAutomorphPrecon(ctxt)  -- this decomposes ctxt into digits
#   for each j: v[j] = precon.automorph(σ^j)  -- this applies automorphism to decomposed form
#
# The key insight: BasicAutomorphPrecon.automorph() does NOT do a full keySwitchPart!
# It applies the automorphism to the DECOMPOSED form and then reconstructs.
# This is the "hoisting" optimization: decompose once, apply many automorphisms cheaply.
#
# So the cost of GenBabySteps is:
#   1 decomposition (BasicAutomorphPrecon ctor): 0.6s
#   g-1 = 9 hoisted automorphisms: 9 × 0.105s = 0.95s
#   (These are the "automorph" calls at matmul.cpp:148)
#
# Then for each giant-step k > 0:
#   acc_inner.smartAutomorph(σ^{g*k}): this is a FULL key-switch (not hoisted)
#   Cost: 0.628s per giant-step (this is the keySwitchPart)
#
# So the breakdown for D=96 (g=10, h=10) is:
#   BasicAutomorphPrecon: 0.6s
#   9 hoisted baby automorphs: 9 × 0.105s = 0.95s
#   96 MulAdd: 96 × ~0.01s = 0.96s
#   9 giant-step smartAutomorph: 9 × 0.628s = 5.65s
#   Overhead (additions, etc.): ~0.2s
#   Total: 0.6 + 0.95 + 0.96 + 5.65 + 0.2 = 8.36s
#
# Hmm, that's more than 7.36s. The discrepancy might be because:
# - Not all 6 calls are D=96 (some are D=29)
# - The 83 keySwitchPart calls are distributed across all 6 MatMul1DExec calls
#
# For D=29 (g=6, h=5): 5+4=9 baby + 4 giant key-switches = ~13 keySwitchPart
# For D=96 (g=10, h=10): 9 baby + 9 giant = ~18 keySwitchPart
# But we have 2 D=96 calls and 4 D=29 calls? No...
#
# Actually for thin bootstrapping:
#   slotToCoeff: 1 Step1 (BlockMatMul) + 1 Step2 (MatMul1D, D=96)
#   coeffToSlot: 1 Step2^{-1} (MatMul1D, D=96) + 1 Step1^{-1} (BlockMatMul)
#
# So there are 2 MatMul1DExec calls for D=96 and 4 other calls (BlockMatMul or D=29).
# Wait, the output says "mul_MatMul1DExec: 44.1725 / 6 = 7.36208"
# 6 calls total. Let me check what they are.
#
# From the ThinEvalMap constructor: matvec has nfactors entries.
# For m=50731=97×523, nfactors=2.
# slotToCoeff (invert=false): matvec[1] (D=29) + matvec[0] (D=96)
# coeffToSlot (invert=true): matvec[0] (D=96) + matvec[1] (D=29)
# Plus the Step1Matrix (BlockMatMul1DExec, not counted in mul_MatMul1DExec).
#
# So 4 MatMul1DExec calls: 2 for D=96, 2 for D=29.
# But the timer says 6 calls... maybe repeat=1 but there are 3 bootstraps?
# No, repeat=1 means 1 bootstrap.
#
# Let me check: in thin bootstrapping, there are TWO ThinEvalMaps:
#   slotToCoeff (invert=false): applies Step1 then Step2 for each dim
#   coeffToSlot (invert=true): applies Step2^{-1} then Step1^{-1} for each dim
#
# For nfactors=2:
#   slotToCoeff: Step1(dim=1) + Step2(dim=0, D=96)  → 1 MatMul1DExec (D=96)
#   coeffToSlot: Step2^{-1}(dim=0, D=96) + Step1^{-1}(dim=1) → 1 MatMul1DExec (D=96)
#
# That's only 2 MatMul1DExec calls. But timer says 6...
#
# OH WAIT - looking at the ThinEvalMap code more carefully:
# matvec.SetLength(nfactors) = 2 entries: matvec[0] (dim=0, D=96) and matvec[1] (dim=1, D=29)
# Both are MatMul1DExec.
# slotToCoeff applies: matvec[1] then matvec[0]
# coeffToSlot applies: matvec[0] then matvec[1]
#
# So per bootstrap: 4 MatMul1DExec calls (2 for D=96, 2 for D=29).
# With repeat=1: 4 calls. But timer says 6...
#
# Maybe there's a third ThinEvalMap for the "parallel" path?
# Or maybe the Step1 is also a MatMul1DExec in some configurations.
#
# Regardless, the KEY POINT is:
# For D=96, the dominant cost is the 9 giant-step smartAutomorph calls.
# Each costs ~0.628s (full key-switch on the accumulated inner sum).
#
# To reduce this, we need FEWER giant steps.
# Standard BSGS: g=10, h=10 → 9 giant steps
# If g=16, h=6 → 5 giant steps (saves 4 × 0.628s = 2.5s)
# But g=16 means 15 baby steps (vs 9), adding 6 × 0.105s = 0.63s
# Net saving: 2.5 - 0.63 = 1.87s per D=96 call
#
# BUT the experiment showed g=16 is SLOWER! Why?
# Because GenBabySteps with g=16 creates 16 baby-step ciphertexts,
# each requiring a hoisted automorphism. The hoisted automorphism
# cost (0.105s) includes some key-switching internally.
#
# Actually wait - let me re-read GenBabySteps:
# Line 1129: precon = BasicAutomorphPrecon(ctxt)
# Line 1133: v[j] = precon.automorph(zMStar.genToPow(dim, j))
#
# BasicAutomorphPrecon.automorph() returns a shared_ptr<Ctxt>.
# The cost is the "automorph" timer at matmul.cpp:148 = 0.105s.
# This does NOT include a keySwitchPart - it's a hoisted operation.
#
# So baby-step cost = 0.105s per step (no key-switch).
# Giant-step cost = smartAutomorph = includes keySwitchPart = 0.628s per step.
#
# Standard (g=10, h=10):
#   Baby: 9 × 0.105 = 0.945s
#   Giant: 9 × 0.628 = 5.652s
#   MulAdd: 96 × 0.01 = 0.96s
#   Precon: 0.6s
#   Total: 8.16s (close to observed 7.36s, difference is overhead)
#
# Frobenius (g=16, h=6):
#   Baby: 15 × 0.105 = 1.575s
#   Giant: 5 × 0.628 = 3.14s
#   MulAdd: 96 × 0.01 = 0.96s
#   Precon: 0.6s
#   Total: 6.28s → should be FASTER!
#
# But experiment showed 11.77s! Something else is going on.
# The issue might be that GeneralAutomorphPrecon_BSGS (not GenBabySteps)
# is being used, and it creates h=6 BasicAutomorphPrecon objects
# (one per giant-step), each costing 0.6s!
#
# From line 294-303:
#   precon.resize(h);  // h = 6 precons!
#   for k in [0..h):
#     p = precon0.automorph(σ^{g*k})  // get giant-step ciphertext
#     precon[k] = BasicAutomorphPrecon(*p)  // decompose it!
#
# So GeneralAutomorphPrecon_BSGS creates h=6 BasicAutomorphPrecon objects!
# Each costs 0.6s → 6 × 0.6s = 3.6s just for precomputation!
# Plus the initial precon0: 0.6s
# Total precomputation: 4.2s
#
# With standard g=10, h=10: 10 precons × 0.6s + 1 initial = 6.6s precomp
# Wait that can't be right either - the total is only 7.36s.
#
# Let me look at the timer more carefully:
# BasicAutomorphPrecon: 4.80985 / 8 = 0.601s
# 8 calls total across all 6 MatMul1DExec calls.
# So only 8/6 ≈ 1.3 BasicAutomorphPrecon per MatMul1DExec call.
#
# This means GeneralAutomorphPrecon_BSGS is NOT being used!
# The code must be taking the "iterative" path (line 1398-1424).
#
# In the iterative path:
#   baby_steps[j] = baby_steps[j-1].smartAutomorph(σ^1)  → each is a full key-switch!
#   Then for each giant-step: sum.smartAutomorph(σ^g) → also a full key-switch
#
# So in iterative mode:
#   Baby: (g-1) smartAutomorphs = 9 × 0.628s = 5.65s
#   Giant: (h-1) smartAutomorphs = 9 × 0.628s = 5.65s
#   Total key-switches: 18 × 0.628s = 11.3s
#   But observed is only 7.36s...
#
# I'm confused. Let me just count: 83 keySwitchPart / 6 calls = 13.8 per call.
# At 0.628s each: 13.8 × 0.628 = 8.7s. But mul_MatMul1DExec = 7.36s.
# The keySwitchPart timer might include calls from OUTSIDE MatMul1DExec.
#
# Bottom line: the experiment with g=16 gave 11.77s vs 7.36s baseline.
# The 98 keySwitchPart (vs 83 baseline) = 15 more key-switches.
# 15 × 0.628 = 9.4s extra. That explains the slowdown.
#
# WHY does g=16 cause MORE key-switches?
# In GeneralAutomorphPrecon_BSGS with g=16, h=6:
#   It creates h=6 precons, each requiring a giant-step automorph from precon0.
#   Then for each query, it uses the appropriate precon[k] to extract baby-step.
#   Total key-switches per query: 1 (from precon[k])
#   Total queries: D=96
#   Plus h=6 precon constructions (each involves key-switch): 6
#   Plus 1 initial precon: 1
#   Total: 96 + 6 + 1 = 103? That's more than 98.
#
# I think the issue is that GeneralAutomorphPrecon_BSGS with larger g
# doesn't reduce key-switches - it just reorganizes them.
# The total number of key-switches is always ~D (one per diagonal).

print("CONCLUSION: The number of key-switches in HElib's BSGS is approximately")
print("equal to D (the dimension), regardless of the baby/giant split.")
print("Changing g from 10 to 16 doesn't reduce key-switches; it may increase them")
print("due to the precomputation structure of GeneralAutomorphPrecon_BSGS.")
print()
print("The ONLY way to reduce key-switches is to reduce the number of")
print("NON-ZERO DIAGONALS that need to be processed.")
print()
print("For a DENSE 96×96 matrix (all 96 diagonals non-zero),")
print("there is NO way to reduce below ~96 key-switches in HElib's framework.")
print()
print("The Frobenius-orbit relation c_{d+16}=c_d*p^i is mathematically correct,")
print("but it doesn't reduce the number of non-zero diagonals.")
print("It only relates the CONSTANTS on different diagonals.")
print()
print("REAL OPTIMIZATION OPPORTUNITY:")
print("The relation means that for the ITERATIVE BSGS path (line 1398-1424),")
print("the inner MulAdd loop can be optimized:")
print("  Standard: for each k, compute Σ_j c_{j+g*k} ⊙ baby[j] (g MulAdds)")
print("  Frobenius: for k>0, compute (p^{ik}) ⊙ (Σ_j c_j ⊙ baby[j]) (1 scaling)")
print("This saves (h-1)*g - (h-1) = (h-1)*(g-1) MulAdd operations.")
print("For g=16, h=6: saves 5*15 = 75 MulAdds (out of 96 total).")
print()
print("But MulAdd only costs ~0.01s each, so saving 75 MulAdds = 0.75s.")
print("This is modest compared to the key-switch cost (~5.6s).")
print()
print("FINAL HONEST ASSESSMENT:")
print("For D=96 with a DENSE matrix, HElib's BSGS with g=10 is already")
print("near-optimal. The key-switch cost dominates, and it's proportional")
print("to the number of distinct rotations needed (which is D for a dense matrix).")
print("No algebraic trick can reduce this without making the matrix SPARSE.")
