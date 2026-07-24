# T2: p=257 order-four measurement in the Ma-NTT power-of-two pipeline, h=7

Param entry added: po2_real_params index i=14, m=2^16=65536, p=257, r=1, c=3, bits=980,
d=ord_65536(257)=256, nslots=128, gens=[-1,5], ords=[2,64], mvec=[65536].
(fatboot.cpp, cloned from the i=6..13 po2_expand template; bits=980 as all i>=6.)

Run: t=-1 newbts=1 newks=1 thick=0 repeat=3, h=7, same-session sequential OFF then ON,
taskset -c 8-25 (cores pinned, matching the po2_grid driver). A=HELIB_EXPLICIT_AUX=16.

## Admissibility (h=7): PASS
bound on I = 6.610101  ->  ceil(I_bound)=7  ->  I_range = 2*7+1 = 15
I_range^2 = 225 < 257  (the recryption.cpp:491 check "p too small to hold the overflow
part" PASSES, slack = 257-225 = 32). aux=16 (coprime to 257). No abort.
security level = 78.39, ord(p)=256, nslots=128.

## Results (extract = digit-extraction step, the part order-4 accelerates)
OFF (generic bounded-support PS, HELIB_EXPLICIT_AUX=16):
  deg of poly for non-power-of-p aux = 223
  extract per repeat: 12.637, 12.557, 12.522 s  -> first 12.64 s, mean 12.57 s
  verdict: ### bts finished, everything ok ###  (decryption CORRECT)

ON (order-4, HELIB_EXPLICIT_AUX=16 HELIB_AUX_ORDER4_EVAL=1):
  TRIGGER: HELIB_AUX_ORDER4_EVAL enabled: deg(P)=223, terms(P)=57, deg(Q)=55, terms(Q)=56
  extract per repeat: 7.947, 7.982, 8.006 s     -> first 7.95 s, mean 7.98 s
  verdict: ### bts finished, everything ok ###  (decryption CORRECT)

## Speedup (extract)
mean: 12.57 / 7.98 = 1.58x    first-repeat: 12.64 / 7.95 = 1.59x
(total bootstrap wall OFF 24.72s -> ON 20.28s; the win is entirely in extract, as
linear1/linear2 are unchanged.)

h=7 worked on the first try; h=6/h=5 fallbacks not needed.

## Files
- fatboot.cpp (i=14 entry added; backup fatboot.cpp.bak_p257)
- build/fatboot (rebuilt)
- results/p257_order4/driver.sh, driver.log
- results/p257_order4/h7_OFF.log, h7_ON.log
