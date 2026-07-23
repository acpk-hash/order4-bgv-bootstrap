/*
 * bench_tower_vs_bsgs.cpp - End-to-end bootstrap timing benchmark
 *
 * Parameters: m=50731, D=96, p=65537
 * Measures: ModRaise, CoeffToSlot (tower), DigitExtract (PS), SlotToCoeff (tower)
 * Reports per-phase timing and operation counts.
 */
#include "params.h"
#include "poly_ring.h"
#include "ciphertext.h"
#include "keys.h"
#include "tower_transform.h"
#include "poly_eval.h"
#include "bootstrap.h"
#include <iostream>
#include <chrono>
#include <iomanip>
#include <set>

using namespace tower_bgv;
using namespace NTL;
using namespace std;

struct Timer {
    chrono::high_resolution_clock::time_point t0;
    void start() { t0 = chrono::high_resolution_clock::now(); }
    double ms() {
        return chrono::duration<double,milli>(
            chrono::high_resolution_clock::now() - t0).count();
    }
    double sec() { return ms() / 1000.0; }
};

int main() {
    cout << "================================================================" << endl;
    cout << "  Tower BGV Bootstrap Benchmark" << endl;
    cout << "  m=50731 (97x523), D=96, p=65537, phi(m)=50112" << endl;
    cout << "================================================================" << endl;
    cout << endl;

    Timer timer, phase_timer;

    // ---- Phase 0: Parameter Setup ----
    cout << "[Phase 0] Parameter setup..." << flush;
    timer.start();
    Params params = Params::create_D96();
    cout << " done" << endl;
    cout << "  m=" << params.m << ", p=" << params.p
         << ", phi(m)=" << params.phim << ", D=" << params.D << endl;
    cout << "  Tower: degrees=[";
    for (size_t i = 0; i < params.tower_degrees.size(); i++)
        cout << params.tower_degrees[i] << (i+1<params.tower_degrees.size()?",":"");
    cout << "], strides=[";
    for (size_t i = 0; i < params.tower_strides.size(); i++)
        cout << params.tower_strides[i] << (i+1<params.tower_strides.size()?",":"");
    cout << "]" << endl;

    // ---- Phase 1: Ring Construction ----
    cout << "\n[Phase 1] Building polynomial ring (Phi_" << params.m << ")..." << flush;
    phase_timer.start();
    PolyRing ring(params);
    double ring_time = phase_timer.sec();
    cout << " " << fixed << setprecision(2) << ring_time << " s" << endl;
    cout << "  deg(Phi_m) = " << deg(ring.phi_m) << endl;

    // ---- Phase 2: Key Generation ----
    cout << "\n[Phase 2] Generating keys..." << endl;
    Keys keys;
    keys.p = params.p;
    // Use a 180-bit modulus (3 x 60-bit primes)
    keys.modulus = ZZ(1);
    for (int i = 0; i < 3; i++) keys.modulus <<= 60;
    cout << "  Q = 2^180 (" << NumBits(keys.modulus) << " bits)" << endl;

    cout << "  Secret key (hwt=64)..." << flush;
    phase_timer.start();
    keys.gen_secret(ring, 64);
    cout << " " << fixed << setprecision(2) << phase_timer.sec() << " s" << endl;

    cout << "  Relinearization key..." << flush;
    phase_timer.start();
    keys.gen_relin_key(ring);
    cout << " " << fixed << setprecision(2) << phase_timer.sec() << " s" << endl;

    cout << "  Tower KS matrices..." << flush;
    phase_timer.start();
    // Generate ALL KS matrices needed for the butterfly diagonals
    // For each stage with radix r and stride s, need rotations at offsets (b-a)*s mod D
    // for a,b in {0..r-1}
    std::set<long> all_offsets;
    for (size_t idx = 0; idx < params.tower_degrees.size(); idx++) {
        long r = params.tower_degrees[idx];
        long s = params.tower_strides[idx];
        for (long a = 0; a < r; a++)
            for (long b = 0; b < r; b++) {
                long d = ((b - a) * s % params.D + params.D) % params.D;
                if (d != 0) all_offsets.insert(d);
            }
    }
    for (long offset : all_offsets) {
        long g_fwd = 1;
        for (long i = 0; i < offset; i++)
            g_fwd = (g_fwd * params.g_ell) % params.m;
        keys.gen_ks_matrix(g_fwd, ring);
    }
    double ks_time = phase_timer.sec();
    cout << " " << fixed << setprecision(2) << ks_time << " s"
         << " (" << keys.ks_matrices.size() << " matrices)" << endl;

    // ---- Phase 3: Tower Transform Construction ----
    cout << "\n[Phase 3] Building tower transform..." << flush;
    phase_timer.start();
    TowerTransform tower(params, ring);
    double tower_build_time = phase_timer.sec();
    cout << " " << fixed << setprecision(2) << tower_build_time << " s" << endl;
    cout << "  Layers: " << tower.layers.size() << endl;
    cout << "  Hoistings: " << tower.count_hoistings() << endl;
    cout << "  Rotations: " << tower.count_rotations() << endl;

    // ---- Phase 4: Encrypt test message ----
    cout << "\n[Phase 4] Encrypting test message..." << flush;
    phase_timer.start();
    ZZX msg;
    SetCoeff(msg, 0, 42);
    SetCoeff(msg, 1, 7);
    SetCoeff(msg, 5, 3);
    Ciphertext ct = encrypt(msg, keys.secret, ring, keys.modulus, params.p);
    cout << " " << fixed << setprecision(2) << phase_timer.sec() << " s" << endl;

    // Verify decryption
    ZZX dec = decrypt(ct, keys.secret, ring);
    cout << "  Decrypt check: coeff[0]=" << coeff(dec, 0)
         << " (expected 42)" << endl;

    // ---- Phase 5: Single rotation timing ----
    cout << "\n[Phase 5] Timing single rotation (key-switch)..." << flush;
    long test_k = keys.ks_matrices.begin()->first;
    phase_timer.start();
    Ciphertext ct_rot = automorph(ct, test_k, keys, ring);
    double single_rot_time = phase_timer.sec();
    cout << " " << fixed << setprecision(2) << single_rot_time << " s" << endl;

    // ---- Phase 6: Single multiplication timing ----
    cout << "\n[Phase 6] Timing single multiplication + relin..." << flush;
    phase_timer.start();
    Ciphertext3 ct3 = mul(ct, ct, ring);
    Ciphertext ct_sq = relinearize(ct3, keys, ring);
    double single_mul_time = phase_timer.sec();
    cout << " " << fixed << setprecision(2) << single_mul_time << " s" << endl;

    // ---- Phase 7: Tower S2C timing ----
    cout << "\n[Phase 7] Tower SlotToCoeff (S2C) - " << tower.layers.size() << " layers..." << endl;
    phase_timer.start();
    Ciphertext ct_s2c = tower.apply(ct, keys);
    double s2c_time = phase_timer.sec();
    cout << "  Tower S2C time: " << fixed << setprecision(2) << s2c_time << " s" << endl;

    // Tower C2S
    cout << "  Tower CoeffToSlot (C2S)..." << flush;
    phase_timer.start();
    Ciphertext ct_c2s = tower.apply_inverse(ct, keys);
    double c2s_time = phase_timer.sec();
    cout << " " << fixed << setprecision(2) << c2s_time << " s" << endl;
    double tower_total = s2c_time + c2s_time;
    cout << "  Tower total (S2C+C2S): " << fixed << setprecision(2) << tower_total << " s" << endl;

    // ---- Phase 8: BSGS linear transform (actual measurement) ----
    cout << "\n[Phase 8] BSGS linear transform (actual, sqrt decomp)..." << endl;
    double bsgs_total = 0;
    {
        // BSGS with baby=10, giant=10 for D=96
        // Generate KS matrices for BSGS rotations: baby steps 1..9, giant steps 10,20,...,90
        long baby = 10, giant = 10;
        cout << "  Generating BSGS KS matrices (baby=" << baby << ", giant=" << giant << ")..." << flush;
        phase_timer.start();
        for (long b = 1; b < baby; b++) {
            long g = 1;
            for (long i = 0; i < b; i++) g = (g * params.g_ell) % params.m;
            if (keys.ks_matrices.find(g) == keys.ks_matrices.end())
                keys.gen_ks_matrix(g, ring);
        }
        for (long gs = 1; gs < giant; gs++) {
            long step = gs * baby;
            long g = 1;
            for (long i = 0; i < step; i++) g = (g * params.g_ell) % params.m;
            if (keys.ks_matrices.find(g) == keys.ks_matrices.end())
                keys.gen_ks_matrix(g, ring);
        }
        cout << " " << fixed << setprecision(2) << phase_timer.sec() << " s"
             << " (total KS: " << keys.ks_matrices.size() << ")" << endl;

        // BSGS linear transform: out = sum_{d=0}^{D-1} diag_d * rot_d(ct)
        // Decompose d = b + g*baby where b in [0,baby), g in [0,giant)
        // out = sum_g rot_{g*baby}( sum_b diag_{b+g*baby} * rot_b(ct) )
        // This requires: baby rotations (hoisted) + giant rotations
        cout << "  Running BSGS S2C..." << flush;
        phase_timer.start();

        // Precompute baby-step rotations (hoisted)
        HoistedDecomp hd_bsgs = hoist(ct, keys, ring);
        std::vector<Ciphertext> baby_rots(baby, Ciphertext(ct.modulus, ct.ptxt_space));
        baby_rots[0] = ct;
        for (long b = 1; b < baby; b++) {
            long g = 1;
            for (long i = 0; i < b; i++) g = (g * params.g_ell) % params.m;
            baby_rots[b] = apply_rotation(ct, hd_bsgs, g, keys, ring);
        }

        // For each giant step, accumulate baby-step contributions then rotate
        Ciphertext bsgs_result(ct.modulus, ct.ptxt_space);
        NTL::clear(bsgs_result.c0);
        NTL::clear(bsgs_result.c1);

        long w_approx = tower.omega_D;  // use same approximate root
        for (long gs = 0; gs < giant; gs++) {
            // Accumulate: inner = sum_b diag_{b+gs*baby} * rot_b(ct)
            Ciphertext inner(ct.modulus, ct.ptxt_space);
            NTL::clear(inner.c0);
            NTL::clear(inner.c1);

            for (long b = 0; b < baby; b++) {
                long d = b + gs * baby;
                if (d >= params.D) break;
                // Build diagonal d: diag_d[i] = w^{i * ((i+d) mod D)}
                std::vector<long> diag_vals(params.D);
                for (long i = 0; i < params.D; i++) {
                    long j = (i + d) % params.D;
                    diag_vals[i] = NTL::PowerMod(w_approx, (long)(((long long)i * j) % params.D), params.p);
                }
                NTL::ZZX diag_poly;
                ring.encode_diagonal(diag_poly, diag_vals, params.p);
                Ciphertext term = mul_const(baby_rots[b], diag_poly, ring);
                inner = add(inner, term, ring);
            }

            // Giant-step rotation
            if (gs == 0) {
                bsgs_result = inner;
            } else {
                long g_giant = 1;
                long step = gs * baby;
                for (long i = 0; i < step; i++) g_giant = (g_giant * params.g_ell) % params.m;
                Ciphertext rotated_inner = automorph(inner, g_giant, keys, ring);
                bsgs_result = add(bsgs_result, rotated_inner, ring);
            }
        }
        double bsgs_s2c_time = phase_timer.sec();
        cout << " " << fixed << setprecision(2) << bsgs_s2c_time << " s" << endl;

        // BSGS C2S (same cost as S2C)
        cout << "  Running BSGS C2S..." << flush;
        phase_timer.start();

        HoistedDecomp hd_bsgs2 = hoist(ct, keys, ring);
        std::vector<Ciphertext> baby_rots2(baby, Ciphertext(ct.modulus, ct.ptxt_space));
        baby_rots2[0] = ct;
        for (long b = 1; b < baby; b++) {
            long g = 1;
            for (long i = 0; i < b; i++) g = (g * params.g_ell) % params.m;
            baby_rots2[b] = apply_rotation(ct, hd_bsgs2, g, keys, ring);
        }

        Ciphertext bsgs_result2(ct.modulus, ct.ptxt_space);
        NTL::clear(bsgs_result2.c0);
        NTL::clear(bsgs_result2.c1);

        long w_inv_approx = NTL::InvMod(w_approx, params.p);
        long D_inv = NTL::InvMod((long)params.D, params.p);
        for (long gs = 0; gs < giant; gs++) {
            Ciphertext inner(ct.modulus, ct.ptxt_space);
            NTL::clear(inner.c0);
            NTL::clear(inner.c1);
            for (long b = 0; b < baby; b++) {
                long d = b + gs * baby;
                if (d >= params.D) break;
                std::vector<long> diag_vals(params.D);
                for (long i = 0; i < params.D; i++) {
                    long j = (i + d) % params.D;
                    long val = NTL::PowerMod(w_inv_approx, (long)(((long long)i * j) % params.D), params.p);
                    diag_vals[i] = (long)(((long long)val * D_inv) % params.p);
                }
                NTL::ZZX diag_poly;
                ring.encode_diagonal(diag_poly, diag_vals, params.p);
                Ciphertext term = mul_const(baby_rots2[b], diag_poly, ring);
                inner = add(inner, term, ring);
            }
            if (gs == 0) {
                bsgs_result2 = inner;
            } else {
                long g_giant = 1;
                long step = gs * baby;
                for (long i = 0; i < step; i++) g_giant = (g_giant * params.g_ell) % params.m;
                Ciphertext rotated_inner = automorph(inner, g_giant, keys, ring);
                bsgs_result2 = add(bsgs_result2, rotated_inner, ring);
            }
        }
        double bsgs_c2s_time = phase_timer.sec();
        cout << "  BSGS C2S: " << fixed << setprecision(2) << bsgs_c2s_time << " s" << endl;

        bsgs_total = bsgs_s2c_time + bsgs_c2s_time;
        cout << "  BSGS total (S2C+C2S): " << fixed << setprecision(2) << bsgs_total << " s" << endl;
    }  // end BSGS block

    // ---- Phase 9: Digit Extraction timing (extrapolated) ----
    // Full PS evaluation of degree-65536 polynomial requires ~512 muls.
    // We measure a few iterations and extrapolate.
    cout << "\n[Phase 9] Digit Extraction timing..." << endl;
    long ps_k = 256;  // ceil(sqrt(65537))
    long ps_m = 256;  // ceil(65537/256)
    cout << "  Paterson-Stockmeyer: degree=" << (params.p-1)
         << ", k=" << ps_k << " baby steps, m=" << ps_m << " giant steps" << endl;
    cout << "  Measuring 3 sample multiplications for extrapolation..." << flush;

    // Time 3 mul+relin operations for accurate per-op estimate
    phase_timer.start();
    Ciphertext ct_tmp = ct;
    for (int i = 0; i < 3; i++) {
        Ciphertext3 tmp3 = mul(ct_tmp, ct, ring);
        ct_tmp = relinearize(tmp3, keys, ring);
    }
    double sample_3mul_time = phase_timer.sec();
    double per_mul_time = sample_3mul_time / 3.0;
    cout << " " << fixed << setprecision(2) << per_mul_time << " s/mul" << endl;

    // Extrapolate: baby steps need k-1=255 muls, giant steps need m-1=255 muls
    // Plus m baby-step evaluations (each ~k scalar muls, cheap)
    double extract_time_est = (ps_k - 1 + ps_m - 1) * per_mul_time;
    cout << "  Estimated digit extraction: " << fixed << setprecision(0)
         << extract_time_est << " s (" << extract_time_est/60.0 << " min)" << endl;
    cout << "  (= " << (ps_k-1) << " baby-step muls + " << (ps_m-1) << " giant-step muls)" << endl;

    // ---- Phase 10: Full Bootstrap (extrapolated) ----
    cout << "\n[Phase 10] Full Bootstrap (measured + extrapolated)..." << endl;
    double total_boot_est = tower_total + extract_time_est;
    cout << "  ModRaise:      ~0 s (trivial)" << endl;
    cout << "  CoeffToSlot:   " << c2s_time << " s (measured)" << endl;
    cout << "  DigitExtract:  " << extract_time_est << " s (extrapolated)" << endl;
    cout << "  SlotToCoeff:   " << s2c_time << " s (measured)" << endl;
    cout << "  TOTAL:         " << total_boot_est << " s ("
         << total_boot_est/60.0 << " min)" << endl;

    // ---- Results Summary ----
    cout << "\n================================================================" << endl;
    cout << "  RESULTS SUMMARY" << endl;
    cout << "================================================================" << endl;
    cout << fixed << setprecision(2);
    cout << "\n  Parameters:" << endl;
    cout << "    Ring dimension N = " << params.phim << endl;
    cout << "    Modulus Q = 2^" << NumBits(keys.modulus) << " (" << NumBits(keys.modulus) << " bits)" << endl;
    cout << "    Plaintext p = " << params.p << endl;
    cout << "    Tower dimension D = " << params.D << endl;
    cout << "    KS digits = " << params.dks << endl;

    cout << "\n  Setup costs:" << endl;
    cout << "    Ring construction:  " << ring_time << " s" << endl;
    cout << "    KS key generation:  " << ks_time << " s" << endl;
    cout << "    Tower build:        " << tower_build_time << " s" << endl;

    cout << "\n  Per-operation costs:" << endl;
    cout << "    Single rotation:    " << single_rot_time << " s" << endl;
    cout << "    Single mul+relin:   " << single_mul_time << " s" << endl;
    cout << "    Avg mul (sampled):  " << per_mul_time << " s" << endl;

    cout << "\n  Bootstrap breakdown:" << endl;
    cout << "    ModRaise:           ~0 s" << endl;
    cout << "    CoeffToSlot (C2S):  " << c2s_time << " s" << endl;
    cout << "    DigitExtract:       " << extract_time_est << " s"
         << " (" << extract_time_est/60.0 << " min) [extrapolated]" << endl;
    cout << "    SlotToCoeff (S2C):  " << s2c_time << " s" << endl;
    cout << "    TOTAL:              " << total_boot_est << " s"
         << " (" << total_boot_est/60.0 << " min)" << endl;

    cout << "\n  Tower vs BSGS comparison (linear transform only):" << endl;
    cout << "    Tower approach (D=96):" << endl;
    cout << "      Layers: " << tower.layers.size() << endl;
    cout << "      Hoistings: " << tower.count_hoistings() << endl;
    cout << "      Rotations: " << tower.count_rotations() << endl;
    cout << "      S2C + C2S time: " << tower_total << " s (measured)" << endl;

    cout << "    BSGS approach (D=96, baby=10, giant=10):" << endl;
    cout << "      Baby steps: 10 (hoisted)" << endl;
    cout << "      Giant steps: 10" << endl;
    cout << "      S2C + C2S time: " << bsgs_total << " s (measured)" << endl;

    double speedup = bsgs_total / tower_total;
    cout << "\n    >>> Linear transform speedup (BSGS/Tower): "
         << fixed << setprecision(2) << speedup << "x <<<" << endl;

    cout << "\n================================================================" << endl;
    cout << "  Benchmark complete." << endl;
    cout << "  Total wall time: " << fixed << setprecision(1) << timer.sec() << " s" << endl;
    cout << "================================================================" << endl;

    return 0;
}
