#!/usr/bin/env python3
"""Plaintext test of the order-four factorisation (paper Sec. 4, Theorems 2-4).

For a prime p, a radix A with A^2 = -1 (mod p) and a bound B, the digit-
extraction polynomial P_A is the interpolant over F_p of

    P_A(h*A + l) = l        for all |h|, |l| <= B,

on the |S_A| = (2B+1)^2 points of the support S_A. This test builds P_A from
that definition alone and checks

  T1  the support has (2B+1)^2 distinct points mod p     (B-injectivity, Thm 3)
  T2  deg P_A = |S_A| - 2                                 (P_A is odd)
  T3  P_A(A X) + A P_A(X) = A X coefficient-wise          (Thm 2)
  T4  nonzero monomials only at k = 1 and k = 3 mod 4     (order-four filter)
  T5  at most (|S_A| - 1)/4 + 1 nonzero terms             (Thm 4; equality unless
      a coefficient vanishes by accident mod p, which happens at p = 2521, B = 13)
  T6  coefficient of X is 1/2, so P_A = X/2 + X^3 Q(X^4)
  T7  P_A and the factored form agree on every point of S_A
  T8  (p = 65537, B = 17 only) P_A and Q equal the polynomial files in 5_results/supplementary/polynomials/,
      coefficient by coefficient
  T9  deg P_A equals the degree HElib printed in every shipped log

It runs on the paper's main point and on every (p, A, B) used in the
ciphertext tables. Pure Python, no dependencies, about 30 s.
"""
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]


def interpolate(xs, ys, p):
    """Coefficients (low to high) of the interpolant mod p, Newton form, O(n^2)."""
    n = len(xs)
    c = list(ys)
    for j in range(1, n):
        for i in range(n - 1, j - 1, -1):
            c[i] = (c[i] - c[i - 1]) * pow(xs[i] - xs[i - j], p - 2, p) % p
    poly = [0] * n
    for k in range(n - 1, -1, -1):          # Horner: poly = poly*(X - x_k) + c_k
        for i in range(n - 1, 0, -1):
            poly[i] = (poly[i - 1] - xs[k] * poly[i]) % p
        poly[0] = (c[k] - xs[k] * poly[0]) % p
    return poly


def evaluate(poly, x, p):
    r = 0
    for a in reversed(poly):
        r = (r * x + a) % p
    return r


def check_point(p, A, B, label, failures, verbose=True):
    notes = NOTES
    def expect(cond, what):
        if not cond:
            failures.append(f"{label}: {what}")
        return cond

    assert A * A % p == p - 1, f"{label}: A^2 != -1"
    pts = {}
    for h in range(-B, B + 1):
        for l in range(-B, B + 1):
            pts[(h * A + l) % p] = l % p
    size = (2 * B + 1) ** 2
    expect(len(pts) == size, f"T1 support has {len(pts)} distinct points, expected {size}")
    xs = list(pts)
    P = interpolate(xs, [pts[x] for x in xs], p)
    while P and P[-1] == 0:
        P.pop()
    deg = len(P) - 1
    expect(deg == size - 2, f"T2 deg P_A = {deg}, expected {size - 2}")
    T = [(P[k] * pow(A, k, p) + A * P[k]) % p for k in range(len(P))]
    expect(T[1] == A % p and all(t == 0 for k, t in enumerate(T) if k != 1),
           "T3 P_A(AX) + A P_A(X) != AX")
    nz = [k for k, a in enumerate(P) if a]
    expect(all(k == 1 or k % 4 == 3 for k in nz), "T4 monomial outside k=1, k=3 mod 4")
    bound = (size - 1) // 4 + 1
    expect(len(nz) <= bound, f"T5 {len(nz)} terms, bound {bound}")
    if len(nz) < bound:
        missing = [k for k in range(3, deg + 1, 4) if not P[k]]
        notes.append(f"{label}: {len(nz)} terms, {bound - len(nz)} below the bound "
                     f"(coefficient of X^{missing[0]} vanishes mod {p})")
    half = pow(2, p - 2, p)
    expect(P[1] == half, f"T6 coefficient of X is {P[1]}, expected 1/2 = {half}")
    Q = [P[4 * j + 3] if 4 * j + 3 < len(P) else 0 for j in range((deg - 3) // 4 + 1)]
    ok7 = all((half * x + pow(x, 3, p) * evaluate(Q, pow(x, 4, p), p)) % p == pts[x] for x in xs)
    expect(ok7, "T7 factored form disagrees with the definition on S_A")
    if verbose:
        print(f"  {label:<26} p={p:<6} A={A:<4} B={B:<3} |S_A|={size:<5} deg P_A={deg:<5} "
              f"terms={len(nz):<4} deg Q={len(Q) - 1:<4} "
              f"{'ok' if not [f for f in failures if f.startswith(label)] else 'FAIL'}")
    return P, Q, deg


NOTES = []


def read_poly(path):
    coeffs = {}
    for line in path.read_text().splitlines():
        if line.strip() and not line.startswith("#"):
            k, a = map(int, line.split())
            coeffs[k] = a
    return coeffs


def main():
    failures = []
    print("Order-four factorisation, plaintext check over F_p")

    # paper's main point, and the two heavier bounds of the degree-scaling sweep
    P, Q, _ = check_point(65537, 256, 17, "m50731 main (h=12)", failures)
    shipped_P = read_poly(ROOT / "5_results" / "supplementary" / "polynomials" / "cleaner_polynomial_p65537_A256_B17.txt")
    shipped_Q = read_poly(ROOT / "5_results" / "supplementary" / "polynomials" / "Q_polynomial_p65537_A256_B17.txt")
    ours_P = {k: a for k, a in enumerate(P) if a}
    ours_Q = {j: a for j, a in enumerate(Q) if a}
    if ours_P != shipped_P or ours_Q != shipped_Q:
        failures.append("T8 interpolant differs from 5_results/supplementary/polynomials/")
    else:
        print(f"  T8 equal to the shipped P_A ({len(shipped_P)} terms) and Q "
              f"({len(shipped_Q)} terms), coefficient by coefficient")

    # every (p, A, B) of the ciphertext tables; B from the degree HElib printed
    logs = sorted((ROOT / "5_results" / "paper_runs").glob("*/*"))
    seen = {}
    deg_re = re.compile(r"deg of poly for non-power-of-p aux is (\d+)")
    p_re = re.compile(r"GTestFatboot: p=(\d+)")
    aux_re = re.compile(r"HELIB_EXPLICIT_AUX=(\d+)")
    for log in logs:
        if not log.is_file():
            continue
        text = log.read_text(errors="replace")
        d, pm, am = deg_re.search(text), p_re.search(text), aux_re.search(text)
        if not (d and pm and am) or "HELIB_AUX_ORDER4_EVAL enabled" not in text:
            continue
        p, A, deg = int(pm.group(1)), int(am.group(1)), int(d.group(1))
        B = (int(round((deg + 2) ** 0.5)) - 1) // 2
        seen.setdefault((p, A, B), []).append((log, deg))
    for (p, A, B), runs in sorted(seen.items()):
        label = f"{runs[0][0].parent.name} p={p}"
        _, _, deg = check_point(p, A, B, label, failures)
        for log, d in runs:
            if d != deg:
                failures.append(f"T9 {log.relative_to(ROOT)} printed deg {d}, interpolant has {deg}")
    print(f"  T9 degree matches HElib's output in {sum(len(r) for r in seen.values())} "
          f"order-four runs over {len(seen)} (p, A, B) points")

    for n in NOTES:
        print("  note: " + n)
    if failures:
        print("\nFAILED")
        for f in failures:
            print("  " + f)
        return 1
    print("\nPASSED")
    return 0


if __name__ == "__main__":
    sys.exit(main())
