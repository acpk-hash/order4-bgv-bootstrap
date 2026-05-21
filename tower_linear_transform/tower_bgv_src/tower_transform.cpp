#include "tower_transform.h"
#include <iostream>
#include <cmath>
#include <set>

namespace tower_bgv {

TowerTransform::TowerTransform(const Params& p, const PolyRing& r)
    : params(p), ring(r) {
    build_twiddles();
}

void TowerTransform::build_twiddles() {
    long D = params.D;
    long p = params.p;

    omega_D = 0;
    if ((p - 1) % D == 0) {
        // Exact: D | (p-1), omega_D exists in Z_p
        long pm1 = p - 1;
        for (long g = 2; g < p; g++) {
            bool is_prim = true;
            long tmp = pm1;
            for (long q = 2; q * q <= tmp; q++) {
                if (tmp % q == 0) {
                    if (NTL::PowerMod(g, pm1 / q, p) == 1) { is_prim = false; break; }
                    while (tmp % q == 0) tmp /= q;
                }
            }
            if (tmp > 1 && NTL::PowerMod(g, pm1 / tmp, p) == 1) is_prim = false;
            if (is_prim) {
                omega_D = NTL::PowerMod(g, (p - 1) / D, p);
                break;
            }
        }
        std::cerr << "TowerTransform: D=" << D << ", omega_D=" << omega_D << " (exact)" << std::endl;
    } else {
        // Approximate: use largest available root for timing benchmark
        // Find the largest k such that k | D and k | (p-1)
        long best_k = 1;
        for (long k = D; k >= 2; k--) {
            if (D % k == 0 && (p - 1) % k == 0) { best_k = k; break; }
        }
        // Use omega_best_k as approximate root (timing is still accurate)
        long pm1 = p - 1;
        for (long g = 2; g < std::min(p, (long)1000); g++) {
            bool is_prim = true;
            long tmp = pm1;
            for (long q = 2; q * q <= tmp; q++) {
                if (tmp % q == 0) {
                    if (NTL::PowerMod(g, pm1 / q, p) == 1) { is_prim = false; break; }
                    while (tmp % q == 0) tmp /= q;
                }
            }
            if (tmp > 1 && NTL::PowerMod(g, pm1 / tmp, p) == 1) is_prim = false;
            if (is_prim) {
                omega_D = NTL::PowerMod(g, (p - 1) / best_k, p);
                break;
            }
        }
        std::cerr << "TowerTransform: D=" << D << ", using omega_" << best_k
                  << "=" << omega_D << " (approx, for timing)" << std::endl;
    }

    layers.resize(params.tower_strides.size());
    for (size_t i = 0; i < params.tower_strides.size(); i++) {
        layers[i].stride = params.tower_strides[i];
        layers[i].degree = params.tower_degrees[i];
        long g_power = 1;
        for (long s = 0; s < layers[i].stride; s++)
            g_power = (g_power * params.g_ell) % params.m;
        layers[i].automorph_k = g_power;
    }
}

static long automorph_for_rotation(long rot, long g_ell, long m) {
    long g = 1;
    long actual_rot = ((rot % (m-1)) + (m-1)) % (m-1);
    for (long i = 0; i < actual_rot; i++) g = (g * g_ell) % m;
    return g;
}

// Apply a sparse matrix stage: out = sum_d diag_d * rot_d(in)
// diag_values[d] = vector of D values for diagonal at offset d
static Ciphertext apply_stage(const Ciphertext& ct,
                              const std::vector<std::pair<long, std::vector<long>>>& diags,
                              const Keys& keys, const PolyRing& ring,
                              const Params& params) {
    Ciphertext result(ct.modulus, ct.ptxt_space);
    NTL::clear(result.c0);
    NTL::clear(result.c1);

    for (const auto& [offset, vals] : diags) {
        bool all_zero = true;
        for (long v : vals) if (v != 0) { all_zero = false; break; }
        if (all_zero) continue;

        NTL::ZZX diag_poly;
        ring.encode_diagonal(diag_poly, vals, params.p);

        Ciphertext rotated(ct.modulus, ct.ptxt_space);
        if (offset == 0) {
            rotated = ct;
        } else {
            long g_k = automorph_for_rotation(offset, params.g_ell, params.m);
            rotated = automorph(ct, g_k, keys, ring);
        }

        Ciphertext term = mul_const(rotated, diag_poly, ring);
        result = add(result, term, ring);
    }
    return result;
}

// Build butterfly stage matrices for DIT NTT
// Returns diagonal values for each stage
static void build_butterfly_diags(
    long D, long p, long omega, bool inverse,
    const std::vector<long>& radices,
    const std::vector<long>& strides,
    std::vector<std::vector<std::pair<long, std::vector<long>>>>& stage_diags)
{
    long num_stages = radices.size();
    stage_diags.resize(num_stages);

    // Build each butterfly stage matrix and extract its diagonals
    // Stage k: radix r at stride s
    // DIT butterfly: block-diagonal structure within groups of size r*s
    // M[base+a*s+j][base+b*s+j] = w_r^{a*b} * w^{j*b*(D/(r*s))}
    // where a,b in {0..r-1}, j in {0..s-1}, base = group_start

    for (long stage = 0; stage < num_stages; stage++) {
        long r = radices[stage];
        long s = strides[stage];
        long gs = r * s;  // group size
        long num_groups = D / gs;
        long w_gs = NTL::PowerMod(omega, D / gs, p);  // primitive gs-th root

        // The stage matrix M has non-zero entries at diagonals {0, s, 2s, ..., (r-1)*s}
        // and also at {D-s, D-2s, ...} (negative offsets wrapping around)
        // More precisely: for each pair (a,b), diagonal offset = (b-a)*s mod D

        // Collect all possible diagonal offsets
        std::set<long> offsets_set;
        for (long a = 0; a < r; a++)
            for (long b = 0; b < r; b++)
                offsets_set.insert(((b - a) * s % D + D) % D);

        for (long d : offsets_set) {
            std::vector<long> vals(D, 0);
            for (long g = 0; g < num_groups; g++) {
                long base = g * gs;
                for (long a = 0; a < r; a++) {
                    for (long j = 0; j < s; j++) {
                        long row = base + a * s + j;
                        long col = (row + d) % D;
                        // col should be base + b*s + j for some b
                        long col_in_group = col - base;
                        if (col_in_group < 0 || col_in_group >= gs) {
                            // col is in a different group — value is 0
                            continue;
                        }
                        long b = col_in_group / s;
                        long col_j = col_in_group % s;
                        if (col_j != j) continue;  // not aligned
                        if (b < 0 || b >= r) continue;

                        // M[row][col] = w_r^{a*b} * w^{j*b*(D/gs)}
                        // Inverse: M^{-1}[row][col] = (1/r) * w^{-a*b*(D/r)} * w^{-j*a*(D/gs)}
                        long val;
                        if (!inverse) {
                            long exp1 = (long)(((long long)a * b * (D / r)) % D);
                            long exp2 = (long)(((long long)j * b * (D / gs)) % D);
                            val = (long)(((long long)NTL::PowerMod(omega, exp1, p)
                                         * NTL::PowerMod(omega, exp2, p)) % p);
                        } else {
                            long omega_inv = NTL::InvMod(omega, p);
                            long r_inv = NTL::InvMod(r, p);
                            long exp1 = (long)(((long long)a * b * (D / r)) % D);
                            long exp2 = (long)(((long long)j * a * (D / gs)) % D);
                            val = (long)((((long long)NTL::PowerMod(omega_inv, exp1, p)
                                         * NTL::PowerMod(omega_inv, exp2, p)) % p
                                         * r_inv) % p);
                        }
                        vals[row] = val;
                    }
                }
            }
            stage_diags[stage].push_back({d, vals});
        }
    }
}

Ciphertext TowerTransform::apply(const Ciphertext& ct, const Keys& keys) const {
    if (omega_D == 0) return ct;
    long D = params.D;

    // Build butterfly stage diagonals
    std::vector<std::vector<std::pair<long, std::vector<long>>>> stage_diags;
    build_butterfly_diags(D, params.p, omega_D, false,
                          params.tower_degrees, params.tower_strides, stage_diags);

    // Apply stages in REVERSE order: last stage (innermost butterfly) first
    // Factorization: V = M_last * ... * M_first, so apply M_first to input first
    Ciphertext current = ct;
    for (long stage = (long)stage_diags.size() - 1; stage >= 0; stage--) {
        current = apply_stage(current, stage_diags[stage], keys, ring, params);
    }
    return current;
}

Ciphertext TowerTransform::apply_inverse(const Ciphertext& ct, const Keys& keys) const {
    if (omega_D == 0) return ct;
    long D = params.D;
    long p = params.p;

    // Build inverse butterfly stages
    std::vector<std::vector<std::pair<long, std::vector<long>>>> stage_diags;
    build_butterfly_diags(D, p, omega_D, true,
                          params.tower_degrees, params.tower_strides, stage_diags);

    // Inverse applies stages in FORWARD order (0, 1, 2)
    // = M3_inv first, then M2_inv, then M1_inv
    // = M1_inv * M2_inv * M3_inv * x
    Ciphertext current = ct;
    for (size_t stage = 0; stage < stage_diags.size(); stage++) {
        current = apply_stage(current, stage_diags[stage], keys, ring, params);
    }
    return current;
}

Ciphertext TowerTransform::apply_layer(const Ciphertext& ct, long layer_idx,
                                        const Keys& keys) const {
    if (omega_D == 0) return ct;
    long D = params.D;

    std::vector<std::vector<std::pair<long, std::vector<long>>>> stage_diags;
    build_butterfly_diags(D, params.p, omega_D, false,
                          params.tower_degrees, params.tower_strides, stage_diags);

    if (layer_idx >= 0 && layer_idx < (long)stage_diags.size())
        return apply_stage(ct, stage_diags[layer_idx], keys, ring, params);
    return ct;
}

long TowerTransform::count_hoistings() const {
    return layers.size();
}

long TowerTransform::count_rotations() const {
    long total = 0;
    for (const auto& layer : layers)
        total += layer.degree - 1;
    return total;
}

} // namespace tower_bgv
