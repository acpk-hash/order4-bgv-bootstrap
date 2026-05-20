#!/usr/bin/env python3
"""Plaintext exact Rader prototype for a primitive-prime EvalMap block.

For the promising HElib set E dimension coming from m_i=97, Step2Matrix has
the form

    y_j = sum_{i=0}^{ell-2} f_i * zeta^{i * g^j},

where ell=97, zeta is a primitive ell-th root of unity, and g generates
(Z/ellZ)^*.  This is exactly a prime-length NTT of the coefficient vector
with the trivial output omitted.  Rader rewrites the non-trivial outputs as a
cyclic convolution of length ell-1, which can then be computed by exact
mixed-radix NTT because ell-1 = 2^5 * 3.

This script verifies equality with the dense Vandermonde transform over a
prime field q where both ell and ell-1 divide q-1.  It also counts modular
adds/muls for the dense and Rader paths.
"""

from __future__ import annotations

import argparse
import math
import random
from dataclasses import dataclass


@dataclass
class Counter:
    adds: int = 0
    muls: int = 0

    @property
    def total(self) -> int:
        return self.adds + self.muls


def is_prime(n: int) -> bool:
    if n < 2:
        return False
    if n % 2 == 0:
        return n == 2
    d = 3
    while d * d <= n:
        if n % d == 0:
            return False
        d += 2
    return True


def find_prime_congruent_one(modulus: int, start_k: int = 1) -> int:
    k = start_k
    while True:
        q = k * modulus + 1
        if is_prime(q):
            return q
        k += 1


def factor(n: int) -> list[int]:
    out: list[int] = []
    d = 2
    while d * d <= n:
        while n % d == 0:
            out.append(d)
            n //= d
        d += 1 if d == 2 else 2
    if n > 1:
        out.append(n)
    return out


def smallest_factor(n: int) -> int:
    if n % 2 == 0:
        return 2
    d = 3
    while d * d <= n:
        if n % d == 0:
            return d
        d += 2
    return n


def primitive_root_prime(p: int) -> int:
    fs = sorted(set(factor(p - 1)))
    for g in range(2, p):
        if all(pow(g, (p - 1) // f, p) != 1 for f in fs):
            return g
    raise RuntimeError(f"no primitive root found for {p}")


def generator_mod_prime(p: int) -> int:
    return primitive_root_prime(p)


def add(a: int, b: int, q: int, ctr: Counter) -> int:
    ctr.adds += 1
    return (a + b) % q


def mul(a: int, b: int, q: int, ctr: Counter) -> int:
    ctr.muls += 1
    return (a * b) % q


def direct_dft(vec: list[int], root: int, q: int, ctr: Counter) -> list[int]:
    n = len(vec)
    out: list[int] = []
    for k in range(n):
        acc = 0
        wk = pow(root, k, q)
        power = 1
        for x in vec:
            acc = add(acc, mul(x, power, q, ctr), q, ctr)
            power = mul(power, wk, q, ctr)
        out.append(acc)
    return out


def mixed_radix_dft(vec: list[int], root: int, q: int, ctr: Counter) -> list[int]:
    n = len(vec)
    if n == 1:
        return vec[:]
    r = smallest_factor(n)
    if r == n:
        return direct_dft(vec, root, q, ctr)
    m = n // r

    inner: list[list[int]] = [[0] * r for _ in range(m)]
    root_r = pow(root, m, q)
    for n2 in range(m):
        sub = [vec[n1 * m + n2] for n1 in range(r)]
        inner[n2] = mixed_radix_dft(sub, root_r, q, ctr)

    for n2 in range(m):
        for k1 in range(r):
            inner[n2][k1] = mul(inner[n2][k1], pow(root, n2 * k1, q), q, ctr)

    out = [0] * n
    root_m = pow(root, r, q)
    for k1 in range(r):
        sub = [inner[n2][k1] for n2 in range(m)]
        vals = mixed_radix_dft(sub, root_m, q, ctr)
        for k2, val in enumerate(vals):
            out[k1 + r * k2] = val
    return out


def cyclic_convolution_ntt(a: list[int], b: list[int], root: int, q: int, ctr: Counter) -> list[int]:
    n = len(a)
    assert len(b) == n
    ahat = mixed_radix_dft(a, root, q, ctr)
    bhat = mixed_radix_dft(b, root, q, ctr)
    chat = [mul(x, y, q, ctr) for x, y in zip(ahat, bhat)]
    inv_root = pow(root, -1, q)
    c = mixed_radix_dft(chat, inv_root, q, ctr)
    inv_n = pow(n, -1, q)
    return [mul(x, inv_n, q, ctr) for x in c]


def dense_primitive_vandermonde(coeffs: list[int], ell: int, zeta: int, gen: int, q: int) -> tuple[list[int], Counter]:
    ctr = Counter()
    out: list[int] = []
    for j in range(ell - 1):
        exponent = pow(gen, j, ell)
        acc = 0
        # Constants are precomputed in HElib, so only count the online
        # constant multiply and add, not the offline power generation.
        for i, coeff in enumerate(coeffs):
            acc = add(acc, mul(coeff, pow(zeta, exponent * i, q), q, ctr), q, ctr)
        out.append(acc)
    return out, ctr


def rader_primitive_vandermonde(coeffs: list[int], ell: int, zeta: int, gen: int, conv_root: int, q: int) -> tuple[list[int], Counter]:
    ctr = Counter()
    n = ell
    m = ell - 1
    x = coeffs[:] + [0]  # full prime-length NTT input; degree ell-1 term is zero

    a = [x[pow(gen, j, n)] for j in range(m)]
    b = [pow(zeta, pow(gen, j, n), q) for j in range(m)]
    b_neg = [b[(-t) % m] for t in range(m)]
    conv = cyclic_convolution_ntt(a, b_neg, conv_root, q, ctr)

    full = [0] * n
    full[0] = sum(x) % q
    for m_idx, val in enumerate(conv):
        k = pow(gen, -m_idx, n)
        full[k] = add(x[0], val, q, ctr)

    return [full[pow(gen, j, n)] for j in range(m)], ctr


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--ell", type=int, default=97)
    parser.add_argument("--seed", type=int, default=7)
    parser.add_argument("--q", type=int)
    args = parser.parse_args()

    ell = args.ell
    assert is_prime(ell), "Rader prototype expects prime ell"
    modulus = math.lcm(ell, ell - 1)
    q = args.q if args.q is not None else find_prime_congruent_one(modulus)
    assert (q - 1) % modulus == 0 and is_prime(q)

    root = primitive_root_prime(q)
    zeta = pow(root, (q - 1) // ell, q)
    conv_root = pow(root, (q - 1) // (ell - 1), q)
    gen = generator_mod_prime(ell)

    random.seed(args.seed)
    coeffs = [random.randrange(q) for _ in range(ell - 1)]

    dense, dense_ctr = dense_primitive_vandermonde(coeffs, ell, zeta, gen, q)
    rader, rader_ctr = rader_primitive_vandermonde(coeffs, ell, zeta, gen, conv_root, q)

    ok = dense == rader
    print(f"ell={ell} q={q} gen_mod_ell={gen} verified={ok}")
    print(f"dense: adds={dense_ctr.adds} muls={dense_ctr.muls} total={dense_ctr.total}")
    print(f"rader+mixed-radix-NTT: adds={rader_ctr.adds} muls={rader_ctr.muls} total={rader_ctr.total}")
    print(f"operation_ratio={rader_ctr.total / dense_ctr.total:.4f}")
    print(f"conv_length={ell - 1} conv_radices={factor(ell - 1)}")
    if not ok:
        for i, (a, b) in enumerate(zip(dense, rader)):
            if a != b:
                raise RuntimeError(f"mismatch at {i}: dense={a}, rader={b}")


if __name__ == "__main__":
    main()
