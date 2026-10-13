# 2_baselines — unmodified third-party sources

Every baseline the paper compares against is included here **as published**,
at the commit used for the paper (the only omission: HElib's
`documentation/Old_Document/`, old design papers and slides that the build does
not use), so that the build needs no
network access and the baselines cannot drift. Nothing in this directory is
edited: `3_experiments/build_all.sh` copies a tree into `build/src/` and applies
patches only to the copy. To check that a tree is unmodified, compare it with
`git archive <commit>` of the upstream repository.

| directory | upstream | commit | licence | role |
|---|---|---|---|---|
| `HElib-3e337a6/` | https://github.com/homenc/HElib | `3e337a66a91a92d49de6a9505340826b0eb71081` | Apache-2.0 (`LICENSE.md`) | base library of every build |
| `BGV-Boot-for-Large-p-83b54d5/` | https://github.com/msh086/BGV-Boot-for-Large-p | `83b54d5` | MIT (`LICENSE`) | Ma et al., EUROCRYPT 2024: baseline of Table 8; its driver is also the base of our driver |
| `bgv-bootstrapping-with-homomorphic-NTT-800bae0/` | https://github.com/msh086/bgv-bootstrapping-with-homomorphic-NTT | `800bae0` | **none published**, see `NOTICE.md` inside | Ma et al., homomorphic NTT for power-of-two rings: baseline of Table 5 |
| `ntl-11.5.1.tar.gz` | https://libntl.org | 11.5.1, SHA-256 `210d06c31306cbc6eaf6814453c56c776d9d8e8df36d74eb306f6a523d1c6a8a` | LGPL-2.1+ | number theory library used by HElib |

## How each baseline is built

`3_experiments/build_all.sh` produces two executables from these trees **without
any of our code**. They form control group G0 in
`3_experiments/compare_baselines.sh`:

| executable | built from | the only change, required by the upstream README |
|---|---|---|
| `build/bin/baseline_ma_largep` | HElib + `BGV-Boot-for-Large-p/bgvboot-final.patch`, driver `fatboot.cpp` | the placeholder `your_path_to_main_dir` in `extractDigits.cpp` is set to the polynomial-cache directory |
| `build/bin/baseline_ma_ntt` | HElib + `bgv-bootstrapping-with-homomorphic-NTT/github_20241003.patch`, driver `fatboot.cpp` | the placeholder `path_of_fatboot.cpp` in `extractDigits.cpp` is set to the polynomial-cache directory |

In both cases `helib_DIR` in the driver's `CMakeLists.txt` is set to the
install location, as the upstream README instructs.

Our own executables start from the same trees and add the diffs in
`1_implementation/patches/` (see `1_implementation/README.md`).

## Baselines not included

* **HElib native bootstrapping** (Table 8, first row): it is the
  same `HElib-3e337a6/` tree with its own Chen–Han digit extraction; at
  `p = 65537` it does not finish within 12 hours, and the paper claims no ratio
  against it.
* **Geelen et al., EUROCRYPT 2023** (Table 8, second row): their
  polynomial generation needs the commercial Magma system at `p = 65537`; the
  paper's row is an extrapolation, marked as such, and is not measured.
* **Okada et al., ASIACRYPT 2023** (Sec. 6.4): their SEAL-based code is not
  redistributed; `5_results/supplementary/okada_comparison/` has our harness
  and logs.
