# Encrypted BSGS-vs-butterfly crossover trend across block dimension D

All runs: HElib fork local/helib_auxradix_opt, single-thread,
enc_crossover (mixed-radix DIF butterfly vs monolithic BSGS MatMul1D of the same
dense DxD DFT matrix), bits=600, same input ciphertext per instance, both paths
verified by decryption (0 slot mismatches). Butterfly = full chain wall; BSGS =
MatMul1DExec.mul() wall (transposed=0; transposed=1 within 1 percent).

Ring construction (mirrors D=576): pick a prime ell so HElib yields a NATIVE hypercube
dimension of size exactly D, and prime p with D | p-1 and D smooth (zeta_D in F_p,
pure mixed-radix DIF, scalar constants). ell=D+1 prime gives dim D directly
(ord_ell(p)=1). D=288 (289 not prime): ell=2017, ord_2017(p)=7, 2016=7*288 -> native
288. Stage data plaintext-verified (chain == row-permuted DFT) by derive_trend.py.

| D    | m      | p        | native dim source         | d_slot | nslots | radices        | butterfly(s) | BSGS mul(s) | bf/BSGS | winner    |
|------|--------|----------|---------------------------|--------|--------|----------------|--------------|-------------|---------|-----------|
| 96   | 50731  | 65537    | Set E dim0 sz96           | 18     | 2784   | [6,4,4]        | 14.9 (b1500) | 14.9        | ~1.0*   | BSGS      |
| 192  | 17177  | 222337   | ell=193, dim0 native 192  | 88     | 192    | [4,4,4,3]      | 3.61         | 2.89        | 1.25    | BSGS      |
| 288  | 14119  | 199873   | ell=2017, dim0 native 288 | 21     | 576    | [4,4,2,3,3]    | 2.79 (2.71)  | 2.70 (2.66) | 1.03    | tie/BSGS  |
| 576  | 51353  | 997057   | ell=577, dim0 native 576  | 44     | 1152   | [4,4,4,3,3]    | 14.32        | 23.51       | 0.61    | butterfly |
| 1152 | 19601  | 3984769  | ell=1153, dim0 native 1152| 16     | 1152   | [4,4,4,3,3,2]  | 5.46         | 9.26        | 0.59    | butterfly |

* D=96 (Set E, m=50731): README reference, butterfly 59.1s vs BSGS 14.9s at bits=1500.

Within-instance butterfly/BSGS ratio is the crossover signal (both paths run
back-to-back on identical ring+ciphertext, machine load cancels): 1.25 (192) ->
1.03 (288) -> 0.61 (576) -> 0.59 (1152). Crossover (ratio=1) just above D~288;
BSGS wins D<=288, butterfly wins D>=576, monotone. Confirms Contribution-2:
O(log D) mixed-radix butterfly overtakes monolithic BSGS as block dim grows,
break-even near D~300. Absolute seconds are NOT comparable across rows (different
d_slot/nslots per row); only per-row ratio is the metric.

Files: derive_trend.py; stages192/288/1152.txt; enc_crossover_trend;
results_trend/run_D{192,288,1152}_*.log (+ D288 _rerun.log).
