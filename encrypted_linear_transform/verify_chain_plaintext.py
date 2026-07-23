#!/usr/bin/env python3
"""Verify stage chain (M2 o M1) == full Step2 dim0 matrix, in F_p[X]/G(X), p=65537, deg G=18.

Sources:
 - stage diags: order4bgv/tower_linear_transform/sage_verification/tower_stage_coefficients.json
 - true matrix: github_order4_cleaner/step2_matrix_dim0_sz96.txt  (HElib dump)
Convention under test: (M v)[i] = sum_k diag_k[i] * v[(i-k) % 96], S = M2 o M1,
with M1 offsets {0,16,...,80} and M2 offsets {0..15}.
"""
import json, random, sys
from pathlib import Path

ROOT = Path("/home/user/experiments")
D = 96
P = 65537

cfg = json.load(open(ROOT/"order4bgv/tower_linear_transform/sage_verification/tower_stage_coefficients.json"))
G = cfg["G"]  # 19 coeffs, degree 18
d = len(G) - 1
assert G[d] == 1 or G[d] != 0
s1_off = cfg["stage1_offsets"]; s2_off = cfg["stage2_offsets"]
s1 = cfg["stage1_diags"]  # [96][6][18]
s2 = cfg["stage2_diags"]  # [96][16][18]

def pmulmod(a, b):
    # multiply two coeff lists mod G mod P
    res = [0]*(2*d-1)
    for i,ai in enumerate(a):
        if ai:
            for j,bj in enumerate(b):
                if bj:
                    res[i+j] = (res[i+j] + ai*bj) % P
    # reduce mod G (monic? check)
    lead = G[d]
    lead_inv = pow(lead, P-2, P)
    for i in range(len(res)-1, d-1, -1):
        c = res[i]
        if c:
            c = c * lead_inv % P
            for j in range(d+1):
                res[i-d+j] = (res[i-d+j] - c*G[j]) % P
    return [x % P for x in res[:d]]

def padd(a, b):
    return [(x+y) % P for x,y in zip(a,b)]

# parse dump
lines = (ROOT/"github_order4_cleaner/step2_matrix_dim0_sz96.txt").read_text().strip().split("\n")
assert int(lines[0]) == D
S = []
for i in range(1, D+1):
    row = []
    for e in lines[i].split("|"):
        e = e.strip()
        assert e.startswith("[") and e.endswith("]")
        cs = [int(c) % P for c in e[1:-1].split()]
        cs = cs + [0]*(d-len(cs))
        row.append(cs[:d])
    assert len(row) == D
    S.append(row)

Sdiag = [[S[i][(i+k) % D] for i in range(D)] for k in range(D)]  # Sdiag[k][i]

random.seed(20260721)
v = [[random.randrange(P) for _ in range(d)] for _ in range(D)]

def apply_diag(diags_by_off, v, sign):
    # diags_by_off: dict offset -> list of 96 field elems
    out = [[0]*d for _ in range(D)]
    for k, dg in diags_by_off.items():
        for i in range(D):
            src = (i - k) % D if sign == "-" else (i + k) % D
            out[i] = padd(out[i], pmulmod(dg[i], v[src]))
    return out

# reference: full matrix as matrix-vector product out[i] = sum_j S[i][j] * v[j]
ref = [[0]*d for _ in range(D)]
for i in range(D):
    acc = [0]*d
    for j in range(D):
        acc = padd(acc, pmulmod(S[i][j], v[j]))
    ref[i] = acc

# also full-diag form with each sign, as sanity for convention
full_minus = apply_diag({k: Sdiag[k] for k in range(D)}, v, "-")
full_plus  = apply_diag({k: Sdiag[k] for k in range(D)}, v, "+")
print("full-diag - == matvec:", full_minus == ref)
print("full-diag + == matvec:", full_plus == ref)

M1 = {off: [s1[i][j] for i in range(D)] for j, off in enumerate(s1_off)}
M2 = {off: [s2[i][j] for i in range(D)] for j, off in enumerate(s2_off)}

for sign in ("-", "+"):
    for order in ("M2oM1", "M1oM2"):
        a, b = (M1, M2) if order == "M2oM1" else (M2, M1)
        w = apply_diag(a, v, sign)
        out = apply_diag(b, w, sign)
        print(f"sign={sign} order={order}: chain == matvec: {out == ref}")
