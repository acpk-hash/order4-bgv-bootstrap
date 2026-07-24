# Measurement note: p = 8101 general-ring speedup (correction)

The originally reported speedup for p = 8101 (m = 23771, h = 12, A = 90) was
3.57x, computed from `p8101_off_r3.log` / `p8101_on_r3.log` as 76.3 s -> 21.4 s
using the first two OFF repeats. A re-examination of those logs shows the three
OFF repeats were 75.1 s, 77.5 s, and 34.5 s: the first two OFF repeats ran while
unrelated compute jobs were loading the machine, and the load ended during the
OFF phase, so the OFF/ON pairing that normally cancels machine load did not
hold for this set.

Clean remeasurement on an idle machine (same binary configuration,
`HELIB_EXPLICIT_AUX=90`, sequential OFF then ON, `repeat=3`, decryption
verified in both runs):

- OFF (generic bounded-support PS): extract = 30.39 / 30.44 / 30.40 s, mean 30.41 s
  (`p8101_off_clean_r3.log`)
- ON (order-four evaluator, trigger line `deg(P)=1087, terms(P)=273, deg(Q)=271`):
  extract = 17.77 / 18.63 / 18.55 s, mean 18.32 s (`p8101_on_clean_r3.log`)

Corrected speedup: 30.41 / 18.32 = **1.66x** (first-repeat pair 30.39 / 17.77 =
1.71x). The corrected value is in line with the other general-ring sets; the
general-ring speedup range should accordingly be read as approximately
1.66-1.90x rather than 1.69-3.57x, and the "exceeds 2x" depth-based explanation
attached to the 3.57x figure does not apply. The capacity accounting of the
clean runs (ON consumes slightly more extract capacity than OFF: 355.8 vs
323.7 bits) matches the paper's appendix noise analysis (+2 levels for the
factored evaluator).

The original loaded logs are retained unmodified in `../moderate_primes/` for
transparency.
