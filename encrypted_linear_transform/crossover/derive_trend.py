#!/usr/bin/env python3
"""Derive mixed-radix DIF butterfly stages for extra crossover D points.
Reuses the exact method of derive_crossover.py (plaintext-verified: stage product
== row-permuted DFT). Writes stages<D>.txt consumed by enc_crossover.
"""
import sys
from sympy import isprime, primitive_root

def derive(D, p, radices, outfile):
    assert isprime(p) and (p-1) % D == 0, (p, D)
    prod = 1
    for r in radices: prod *= r
    assert prod == D, (prod, D)
    g = primitive_root(p)
    zeta = pow(g, (p-1)//D, p)
    assert pow(zeta, D, p) == 1
    for q in set(radices):
        assert pow(zeta, D//q, p) != 1
    stages = []
    L = D
    for r in radices:
        m2 = L // r
        M = {}
        for b in range(0, D, L):
            for j in range(m2):
                for i in range(r):
                    row = b + j + i*m2
                    tw = pow(zeta, (D//L)*i*j, p)
                    for i2 in range(r):
                        col = b + j + i2*m2
                        M[(row, col)] = pow(zeta, (D//r)*i*i2, p) * tw % p
        stages.append(M)
        L = m2
    assert L == 1
    offs_all = []; nrot = 0
    for si, M in enumerate(stages):
        offs = sorted(set((c - r0) % D for (r0, c) in M))
        offs_all.append(offs)
        nrot += sum(1 for o in offs if o != 0)
        print(f"stage {si+1} (radix {radices[si]}): {len(offs)} diagonals, offsets {offs}")
    print(f"total nontrivial rotations: {nrot}, const-mult levels: {len(stages)}")
    cols = []
    for j in range(D):
        v = [0]*D; v[j] = 1
        for M in stages:
            out = [0]*D
            for (r0, c), val in M.items():
                if v[c]: out[r0] = (out[r0] + val * v[c]) % p
            v = out
        cols.append(v)
    Wrow_index = {}
    for i in range(D):
        Wrow_index[tuple(pow(zeta, i*j, p) for j in range(D))] = i
    perm = []
    for t in range(D):
        row = tuple(cols[j][t] for j in range(D))
        idx = Wrow_index.get(row)
        assert idx is not None, f"chain row {t} is not a DFT row"
        perm.append(idx)
    assert sorted(perm) == list(range(D))
    print("VERIFIED: chain product == row-permuted DFT, perm sample", perm[:8])
    with open(outfile, "w") as f:
        f.write(f"{p} {D} {zeta}\n")
        f.write(f"{len(stages)}\n")
        for si, M in enumerate(stages):
            offs = offs_all[si]
            f.write(f"{len(offs)}\n")
            for o in offs:
                dv = [0]*D
                for (r0, c), val in M.items():
                    if (c - r0) % D == o: dv[r0] = val
                f.write(f"{o}\n")
                f.write(" ".join(map(str, dv)) + "\n")
        f.write(" ".join(map(str, perm)) + "\n")
    print("wrote", outfile)

BASE = "/home/user/experiments/order4bgv/tower_linear_transform/encrypted_prototype/crossover/"
if __name__ == "__main__":
    which = sys.argv[1] if len(sys.argv)>1 else "all"
    if which in ("192","all"):
        print("### D=192 (m=17177=193*89, p=222337) ###")
        derive(192, 222337, [4,4,4,3], BASE+"stages192.txt")
        print()
    if which in ("288","all"):
        print("### D=288 (m=14119=2017*7, p=199873) ###")
        derive(288, 199873, [4,4,2,3,3], BASE+"stages288.txt")
        print()
