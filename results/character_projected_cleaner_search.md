# Character-Projected Cleaner Search

This note records the broader radix search motivated by the matrix-arithmetic
paper and by the finite-field Gaussian-unit interpretation of the implemented
`A=256` cleaner.

## Command

```bash
python3 scripts/character_projected_cleaner_search.py \
  --p 65537 --B 17 --max-order 64 --max-decomposition 16 \
  --aux 35 256 65281 \
  --output results/character_projected_p65537_B17.json --json
```

## General Parameter Rule

For exact cleaner correctness, `p` can be any sufficiently large prime.  A
simple universal choice is:

```text
A = 2B + 1,     p > 4B(B+1).
```

This makes the support map `(h,lo) -> hA+lo mod p` injective on
`|h|,|lo| <= B`, so the exact cleaner exists uniquely as the degree-`<|S_A|`
interpolation polynomial.

For the order-four sparse cleaner used by the optimized HElib experiment,
`p` is not arbitrary.  A clean sufficient theorem is:

```text
p prime, p = 1 mod 4, p > 8B^2, choose A with A^2 = -1 mod p.
```

Then the support is automatically injective, multiplication by `A` rotates
`(h,lo)` to `(lo,-h)`, and the exact cleaner satisfies:

```text
P_A(X) = (1/2) X + X^3 Q(X^4).
```

For `B=17`, the two lower bounds are:

```text
generic correctness: p > 1224
order-four theorem: p > 2312 and p = 1 mod 4
```

The target `p=65537` satisfies both.

## Search Result

The search considered 65 candidate radices, with 44 valid support-injective
candidates and 21 rejected because of support collisions.

| radix `A` | field relation | degree | terms | best decomposition | estimated mults |
|---:|---|---:|---:|---|---:|
| 256 | `A^2=-1` | 1223 | 307 | mod 4: `X/2 + X^3 Q(X^4)` | 41 |
| 65281 | `A^2=-1`, centered `-256` | 1223 | 307 | mod 4 | 41 |
| 35 | generic Ma radix | 1223 | 612 | mod 2 odd-polynomial form | 55 |
| 64 | low-order candidate | 1223 | 612 | mod 2 odd-polynomial form | 55 |

The field-theoretic roots `256` and `65281=-256 mod 65537` are equivalent for
the polynomial identity.  The implementation uses the small positive integer
representative `256`, because the HElib code also uses the auxiliary radix as
an integer scaling parameter.

## Conclusion

The broader search did not find a better implemented candidate than `A=256`.
It did strengthen the paper claim: the improvement is not a one-off SAT
artifact, but a character-projected exact-cleaner construction.  The correct
claim boundary is:

```text
global search/optimization inside the fixed bounded-support exact-cleaner
problem, not global optimality over all bootstrapping algorithms.
```
