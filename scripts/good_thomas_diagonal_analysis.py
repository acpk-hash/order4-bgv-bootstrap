#!/usr/bin/env python3
"""
Compute the exact diagonal structure of the two Good-Thomas stages
after folding permutations into the matrices.

This determines the exact automorphism count for the C++ implementation.
"""
import numpy as np


def power_mod(base, exp, mod):
    result = 1
    base = base % mod
    while exp > 0:
        if exp % 2 == 1:
            result = (result * base) % mod
        exp //= 2
        base = (base * base) % mod
    return result


def find_generator(m):
    phi = m - 1
    factors = []
    n = phi
    for p_f in range(2, int(n**0.5) + 1):
        if n % p_f == 0:
            factors.append(p_f)
            while n % p_f == 0:
                n //= p_f
    if n > 1:
        factors.append(n)
    for g in range(2, m):
        if all(power_mod(g, phi // f, m) != 1 for f in factors):
            return g
    return None


def main():
    p = 65537
    ell = 97
    D = 96
    N1 = 32  # 32-point sub-DFT
    N2 = 3   # 3-point sub-DFT
    cofactor = 523

    g_ell = find_generator(ell)
    cof_mod_ell = cofactor % ell

    # Build the full Step2Matrix (scalar version over F_ell for structure analysis)
    # A[i,j] = (cof * g^j)^i mod ell = cof^i * g^{ij} mod ell
    omega = g_ell  # primitive 96th root mod 97

    A_full = np.zeros((D, D), dtype=np.int64)
    for i in range(D):
        for j in range(D):
            A_full[i][j] = power_mod(cof_mod_ell * power_mod(g_ell, j, ell), i, ell)

    # Good-Thomas permutations
    inv_N2_N1 = pow(N2, -1, N1)  # inv(3,32) = 11
    inv_N1_N2 = pow(N1, -1, N2)  # inv(32,3) = 2

    def crt_reconstruct(k1, k2):
        return (k1 * N2 * inv_N2_N1 + k2 * N1 * inv_N1_N2) % D

    def ruritanian(n1, n2):
        return (n1 * N2 + n2 * N1) % D

    # Build permutation matrices
    P_in = np.zeros((D, D), dtype=np.int64)   # maps linear col to (n1,n2) ordering
    P_out = np.zeros((D, D), dtype=np.int64)  # maps (k1,k2) ordering to linear row

    for n1 in range(N1):
        for n2 in range(N2):
            linear_col = ruritanian(n1, n2)
            tensor_col = n1 * N2 + n2  # row-major (n1, n2)
            P_in[tensor_col][linear_col] = 1

    for k1 in range(N1):
        for k2 in range(N2):
            linear_row = crt_reconstruct(k1, k2)
            tensor_row = k1 * N2 + k2  # row-major (k1, k2)
            P_out[linear_row][tensor_row] = 1

    # Build tensor product (V_32 ⊗ V_3) in tensor indexing
    omega_32 = power_mod(omega, N2, ell)  # ω^3, order 32
    omega_3 = power_mod(omega, N1, ell)   # ω^32, order 3

    T = np.zeros((D, D), dtype=np.int64)
    for k1 in range(N1):
        for k2 in range(N2):
            for n1 in range(N1):
                for n2 in range(N2):
                    row = k1 * N2 + k2
                    col = n1 * N2 + n2
                    T[row][col] = (power_mod(omega_32, k1*n1, ell) *
                                   power_mod(omega_3, k2*n2, ell)) % ell

    # Build diagonal scaling
    D_scale = np.zeros((D, D), dtype=np.int64)
    for i in range(D):
        D_scale[i][i] = power_mod(cof_mod_ell, i, ell)

    # Full factored matrix: D_scale * P_out * T * P_in
    M_factored = D_scale @ P_out @ T @ P_in % ell

    # Verify
    errors = np.sum(M_factored != A_full)
    print(f"Verification: factored vs direct Step2Matrix errors = {errors}")
    assert errors == 0, "Factorization mismatch!"

    # Now build the two-stage decomposition
    # Strategy: fold permutations so both stages operate in LINEAR indexing (0..95)
    #
    # Option A: Stage1 does 3-point, Stage2 does 32-point
    #   M1 = P_mid * (I_32 ⊗ V_3) * P_in   (3-point DFT + input perm)
    #   M2 = D_scale * P_out * (V_32 ⊗ I_3) * P_mid^{-1}  (32-point DFT + output perm + scaling)
    #   where P_mid maps tensor (n1, k2) ordering to linear
    #
    # Simplest: let P_mid = identity in tensor ordering, i.e., intermediate uses tensor index directly
    # Then M1 maps linear → tensor, M2 maps tensor → linear
    #
    # But HElib needs both stages to operate on the SAME linear indexing.
    # Solution: choose P_mid such that intermediate is also in linear indexing.
    #
    # Best choice: P_mid = Ruritanian on (n1, k2) → linear
    # i.e., intermediate linear index = ruritanian(n1, k2) = n1*3 + k2*32 mod 96

    # Actually, the simplest correct approach:
    # Just compute M1 and M2 as full 96×96 matrices in linear indexing,
    # where M2 * M1 = A_full, and M1 is sparse, M2 is moderately sparse.

    # Stage 1: 3-point DFT (transforms n2 → k2, keeps n1 fixed)
    # In tensor indexing: S1_tensor[(n1,k2), (n1',n2)] = δ(n1,n1') * V_3[k2,n2]
    S1_tensor = np.zeros((D, D), dtype=np.int64)
    for n1 in range(N1):
        for k2 in range(N2):
            for n2 in range(N2):
                row = n1 * N2 + k2
                col = n1 * N2 + n2
                S1_tensor[row][col] = power_mod(omega_3, k2 * n2, ell)

    # Stage 2: 32-point DFT (transforms n1 → k1, keeps k2 fixed) + scaling
    # In tensor indexing: S2_tensor[(k1,k2), (n1,k2')] = δ(k2,k2') * V_32[k1,n1]
    S2_tensor = np.zeros((D, D), dtype=np.int64)
    for k1 in range(N1):
        for k2 in range(N2):
            for n1 in range(N1):
                row = k1 * N2 + k2
                col = n1 * N2 + k2  # k2' must equal k2
                S2_tensor[row][col] = power_mod(omega_32, k1 * n1, ell)

    # Verify in tensor space: D_scale * P_out * S2_tensor * S1_tensor * P_in = A_full
    check = D_scale @ P_out @ S2_tensor @ S1_tensor @ P_in % ell
    assert np.sum(check != A_full) == 0, "Tensor stage decomposition failed!"
    print("Tensor stage decomposition verified ✓")

    # Convert to linear indexing:
    # M1 (linear→linear) = P_out * S1_tensor * P_in  ... no, need intermediate
    # Let's use: intermediate = tensor ordering (row-major n1*3+k2)
    # M1_linear[tensor_row, linear_col] = S1_tensor[tensor_row, tensor_col] * P_in[tensor_col, linear_col]
    # = (S1_tensor * P_in)[tensor_row, linear_col]
    M1_linear = S1_tensor @ P_in % ell  # maps linear input → tensor intermediate

    # M2_linear[linear_row, tensor_col] = (D_scale * P_out * S2_tensor)[linear_row, tensor_col]
    M2_linear = D_scale @ P_out @ S2_tensor % ell  # maps tensor intermediate → linear output

    # Verify: M2_linear * M1_linear = A_full
    check2 = M2_linear @ M1_linear % ell
    assert np.sum(check2 != A_full) == 0, "Linear stage decomposition failed!"
    print("Linear stage decomposition verified ✓")

    # But wait - M1_linear maps from linear indexing to tensor indexing,
    # and M2_linear maps from tensor indexing to linear indexing.
    # For HElib, BOTH stages must use the same indexing (linear).
    # We need an intermediate permutation.
    #
    # Solution: define P_mid that maps tensor → linear (any bijection works)
    # Then: M1_hlib = P_mid * S1_tensor * P_in  (linear → linear)
    #        M2_hlib = D_scale * P_out * S2_tensor * P_mid^{-1}  (linear → linear)
    #
    # Choose P_mid = P_out (reuse CRT reconstruction for intermediate)
    # Then: M1_hlib = P_out * S1_tensor * P_in
    #        M2_hlib = D_scale * P_out * S2_tensor * P_out^{-1}

    P_out_inv = np.linalg.inv(P_out.astype(float)).round().astype(np.int64) % ell

    # Actually for permutation matrices, inverse = transpose
    P_out_inv = P_out.T

    M1_hlib = P_out @ S1_tensor @ P_in % ell
    M2_hlib = D_scale @ (P_out @ S2_tensor @ P_out_inv % ell) % ell

    check3 = M2_hlib @ M1_hlib % ell
    assert np.sum(check3 != A_full) == 0, "HElib-compatible decomposition failed!"
    print("HElib-compatible (linear→linear) decomposition verified ✓")

    # Analyze diagonal structure of M1_hlib and M2_hlib
    print("\n=== Diagonal Analysis ===")

    def get_nonzero_diagonals(M, D, mod):
        """Get non-zero diagonal offsets of matrix M (mod arithmetic)."""
        nz_diags = []
        for d in range(D):
            is_nonzero = False
            for i in range(D):
                j = (i + d) % D
                if M[i][j] % mod != 0:
                    is_nonzero = True
                    break
            if is_nonzero:
                nz_diags.append(d)
        return nz_diags

    diags_M1 = get_nonzero_diagonals(M1_hlib, D, ell)
    diags_M2 = get_nonzero_diagonals(M2_hlib, D, ell)

    print(f"\nStage 1 (M1_hlib): {len(diags_M1)} non-zero diagonals")
    print(f"  Offsets: {diags_M1}")

    print(f"\nStage 2 (M2_hlib): {len(diags_M2)} non-zero diagonals")
    if len(diags_M2) <= 40:
        print(f"  Offsets: {diags_M2}")

    # BSGS estimates
    import math
    def bsgs_auts(num_diags):
        if num_diags <= 1:
            return 0
        baby = int(math.ceil(math.sqrt(num_diags)))
        giant = int(math.ceil(num_diags / baby))
        return baby + giant - 1

    auts_M1 = bsgs_auts(len(diags_M1))
    auts_M2 = bsgs_auts(len(diags_M2))

    print(f"\nAutomorphism estimates:")
    print(f"  Stage 1: {auts_M1} (for {len(diags_M1)} diagonals)")
    print(f"  Stage 2: {auts_M2} (for {len(diags_M2)} diagonals)")
    print(f"  Total: {auts_M1 + auts_M2}")
    print(f"  Current BSGS (D=96): 18")
    print(f"  Reduction: {18 - (auts_M1+auts_M2)} fewer ({(18-(auts_M1+auts_M2))/18*100:.1f}%)")

    # Try alternative P_mid choices to minimize total diagonals
    print("\n=== Trying alternative intermediate permutations ===")

    # Try P_mid = Ruritanian (n1,k2) → n1*3 + k2*32 mod 96
    def build_perm_matrix(perm_func, D, N1, N2):
        P = np.zeros((D, D), dtype=np.int64)
        for a in range(N1):
            for b in range(N2):
                tensor_idx = a * N2 + b
                linear_idx = perm_func(a, b)
                P[linear_idx][tensor_idx] = 1
        return P

    # Ruritanian as P_mid
    P_mid_rur = build_perm_matrix(ruritanian, D, N1, N2)
    P_mid_rur_inv = P_mid_rur.T

    M1_v2 = P_mid_rur @ S1_tensor @ P_in % ell
    M2_v2 = D_scale @ P_out @ S2_tensor @ P_mid_rur_inv % ell

    check_v2 = M2_v2 @ M1_v2 % ell
    assert np.sum(check_v2 != A_full) == 0

    diags_M1_v2 = get_nonzero_diagonals(M1_v2, D, ell)
    diags_M2_v2 = get_nonzero_diagonals(M2_v2, D, ell)
    auts_v2 = bsgs_auts(len(diags_M1_v2)) + bsgs_auts(len(diags_M2_v2))
    print(f"P_mid=Ruritanian: Stage1={len(diags_M1_v2)} diags, Stage2={len(diags_M2_v2)} diags, total auts={auts_v2}")

    # Try P_mid = identity (tensor ordering IS linear ordering)
    P_mid_id = np.eye(D, dtype=np.int64)
    M1_v3 = S1_tensor @ P_in % ell
    M2_v3 = D_scale @ P_out @ S2_tensor % ell

    check_v3 = M2_v3 @ M1_v3 % ell
    assert np.sum(check_v3 != A_full) == 0

    diags_M1_v3 = get_nonzero_diagonals(M1_v3, D, ell)
    diags_M2_v3 = get_nonzero_diagonals(M2_v3, D, ell)
    auts_v3 = bsgs_auts(len(diags_M1_v3)) + bsgs_auts(len(diags_M2_v3))
    print(f"P_mid=Identity:   Stage1={len(diags_M1_v3)} diags, Stage2={len(diags_M2_v3)} diags, total auts={auts_v3}")

    # Try swapping: Stage1 does 32-point, Stage2 does 3-point
    # S1'_tensor[(k1,n2), (n1,n2')] = δ(n2,n2') * V_32[k1,n1]
    S1_swap = np.zeros((D, D), dtype=np.int64)
    for k1 in range(N1):
        for n2 in range(N2):
            for n1 in range(N1):
                row = k1 * N2 + n2
                col = n1 * N2 + n2
                S1_swap[row][col] = power_mod(omega_32, k1 * n1, ell)

    # S2'_tensor[(k1,k2), (k1',n2)] = δ(k1,k1') * V_3[k2,n2]
    S2_swap = np.zeros((D, D), dtype=np.int64)
    for k1 in range(N1):
        for k2 in range(N2):
            for n2 in range(N2):
                row = k1 * N2 + k2
                col = k1 * N2 + n2
                S2_swap[row][col] = power_mod(omega_3, k2 * n2, ell)

    # Verify
    check_swap = D_scale @ P_out @ S2_swap @ S1_swap @ P_in % ell
    assert np.sum(check_swap != A_full) == 0, "Swapped order failed!"

    M1_swap = S1_swap @ P_in % ell
    M2_swap = D_scale @ P_out @ S2_swap % ell
    diags_M1_swap = get_nonzero_diagonals(M1_swap, D, ell)
    diags_M2_swap = get_nonzero_diagonals(M2_swap, D, ell)
    auts_swap = bsgs_auts(len(diags_M1_swap)) + bsgs_auts(len(diags_M2_swap))
    print(f"Swapped(32→3):   Stage1={len(diags_M1_swap)} diags, Stage2={len(diags_M2_swap)} diags, total auts={auts_swap}")

    # Find the best configuration
    configs = [
        ("P_mid=CRT", len(diags_M1), len(diags_M2), auts_M1 + auts_M2),
        ("P_mid=Rur", len(diags_M1_v2), len(diags_M2_v2), auts_v2),
        ("P_mid=Id", len(diags_M1_v3), len(diags_M2_v3), auts_v3),
        ("Swapped", len(diags_M1_swap), len(diags_M2_swap), auts_swap),
    ]
    configs.sort(key=lambda x: x[3])
    print(f"\n=== Best configuration: {configs[0][0]} with {configs[0][3]} total auts ===")
    print(f"  Stage 1: {configs[0][1]} diagonals")
    print(f"  Stage 2: {configs[0][2]} diagonals")

    # Print the best stage diagonal offsets for C++ implementation
    best = configs[0][0]
    if best == "P_mid=Id":
        print(f"\n  Stage 1 diagonal offsets: {diags_M1_v3}")
        print(f"  Stage 2 diagonal offsets: {diags_M2_v3}")
    elif best == "Swapped":
        print(f"\n  Stage 1 diagonal offsets: {diags_M1_swap}")
        print(f"  Stage 2 diagonal offsets: {diags_M2_swap}")


if __name__ == "__main__":
    main()
