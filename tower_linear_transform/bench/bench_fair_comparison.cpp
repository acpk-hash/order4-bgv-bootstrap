/*
 * bench_fair_comparison.cpp
 * Fair comparison: HElib-style (BSGS) bootstrap vs Tower-optimized bootstrap
 * Same parameters, same arithmetic, same digit extraction — only linear transform differs.
 */
#include "params.h"
#include "poly_ring.h"
#include "ciphertext.h"
#include "keys.h"
#include "tower_transform.h"
#include <iostream>
#include <chrono>
#include <iomanip>
#include <set>

using namespace tower_bgv;
using namespace NTL;
using namespace std;

typedef std::set<long> LongSet;

struct Timer {
    chrono::high_resolution_clock::time_point t0;
    void start() { t0 = chrono::high_resolution_clock::now(); }
    double sec() { return chrono::duration<double>(chrono::high_resolution_clock::now()-t0).count(); }
};

int main() {
    cout << "================================================================" << endl;
    cout << "  Fair Bootstrap Comparison: HElib-style vs Tower" << endl;
    cout << "  Same arithmetic, same parameters, same digit extraction" << endl;
    cout << "  Only the linear transform algorithm differs" << endl;
    cout << "================================================================\n" << endl;

    Timer timer, phase;
    timer.start();

    // ---- Setup ----
    Params params = Params::create_D96();
    cout << "[Setup] m=" << params.m << ", p=" << params.p
         << ", N=" << params.phim << ", D=" << params.D << endl;

    cout << "  Building ring..." << flush;
    phase.start();
    PolyRing ring(params);
    cout << " " << fixed << setprecision(1) << phase.sec() << "s" << endl;

    Keys keys;
    keys.p = params.p;
    keys.modulus = ZZ(1) << 180;
    keys.gen_secret(ring, 64);

    cout << "  Generating relin key..." << flush;
    phase.start();
    keys.gen_relin_key(ring);
    cout << " " << fixed << setprecision(1) << phase.sec() << "s" << endl;

    // Generate ALL KS matrices needed for both approaches
    cout << "  Generating KS matrices..." << flush;
    phase.start();

    // Tower needs: all butterfly diagonal offsets
    LongSet tower_offsets;
    for (size_t idx = 0; idx < params.tower_degrees.size(); idx++) {
        long r = params.tower_degrees[idx];
        long s = params.tower_strides[idx];
        for (long a = 0; a < r; a++)
            for (long b = 0; b < r; b++) {
                long d = ((b-a)*s % params.D + params.D) % params.D;
                if (d != 0) tower_offsets.insert(d);
            }
    }

    // BSGS needs: baby steps 1..9, giant steps 10,20,...,90
    long baby = 10, giant = 10;
    LongSet bsgs_offsets;
    for (long b = 1; b < baby; b++) bsgs_offsets.insert(b);
    for (long g = 1; g < giant; g++) bsgs_offsets.insert(g * baby);

    // Union of all offsets
    LongSet all_offsets;
    all_offsets.insert(tower_offsets.begin(), tower_offsets.end());
    all_offsets.insert(bsgs_offsets.begin(), bsgs_offsets.end());

    for (long offset : all_offsets) {
        long g = 1;
        for (long i = 0; i < offset; i++) g = (g * params.g_ell) % params.m;
        keys.gen_ks_matrix(g, ring);
    }
    cout << " " << fixed << setprecision(1) << phase.sec() << "s"
         << " (" << keys.ks_matrices.size() << " matrices)" << endl;

    // Encrypt test message
    ZZX msg; SetCoeff(msg, 0, 42);
    Ciphertext ct = encrypt(msg, keys.secret, ring, keys.modulus, params.p);

    // ---- Measure digit extraction cost (same for both) ----
    cout << "\n[Digit Extraction] Measuring per-multiplication cost..." << flush;
    phase.start();
    Ciphertext ct_tmp = ct;
    const int NUM_SAMPLE_MULS = 5;
    for (int i = 0; i < NUM_SAMPLE_MULS; i++) {
        Ciphertext3 tmp3 = mul(ct_tmp, ct, ring);
        ct_tmp = relinearize(tmp3, keys, ring);
    }
    double sample_time = phase.sec();
    double per_mul = sample_time / NUM_SAMPLE_MULS;
    long ps_muls = 510;  // PS for degree 65536: k=256 baby + 255 giant = 510 muls
    double de_time = ps_muls * per_mul;
    cout << " " << fixed << setprecision(2) << per_mul << " s/mul" << endl;
    cout << "  Paterson-Stockmeyer: " << ps_muls << " muls × "
         << per_mul << " s = " << fixed << setprecision(0) << de_time << " s ("
         << de_time/60 << " min)" << endl;

    // ---- BASELINE: HElib-style BSGS Bootstrap ----
    cout << "\n================================================================" << endl;
    cout << "  BASELINE: HElib-style Bootstrap (BSGS linear transform)" << endl;
    cout << "================================================================" << endl;

    // BSGS S2C
    cout << "\n  [1/3] SlotToCoeff (BSGS, baby=" << baby << " giant=" << giant << ")..." << flush;
    phase.start();
    {
        HoistedDecomp hd = hoist(ct, keys, ring);
        vector<Ciphertext> baby_rots(baby, Ciphertext(ct.modulus, ct.ptxt_space));
        baby_rots[0] = ct;
        for (long b = 1; b < baby; b++) {
            long g = 1;
            for (long i = 0; i < b; i++) g = (g * params.g_ell) % params.m;
            baby_rots[b] = apply_rotation(ct, hd, g, keys, ring);
        }

        Ciphertext result(ct.modulus, ct.ptxt_space);
        NTL::clear(result.c0); NTL::clear(result.c1);
        long w = 65529;  // omega_32 (approximate twiddle for timing)
        for (long gs = 0; gs < giant; gs++) {
            Ciphertext inner(ct.modulus, ct.ptxt_space);
            NTL::clear(inner.c0); NTL::clear(inner.c1);
            for (long b = 0; b < baby; b++) {
                long d = b + gs * baby;
                if (d >= params.D) break;
                vector<long> diag_vals(params.D);
                for (long i = 0; i < params.D; i++) {
                    long j = (i + d) % params.D;
                    diag_vals[i] = PowerMod(w, (long)(((long long)i*j) % params.D), params.p);
                }
                ZZX diag_poly; ring.encode_diagonal(diag_poly, diag_vals, params.p);
                Ciphertext term = mul_const(baby_rots[b], diag_poly, ring);
                inner = add(inner, term, ring);
            }
            if (gs == 0) { result = inner; }
            else {
                long g = 1; long step = gs*baby;
                for (long i = 0; i < step; i++) g = (g * params.g_ell) % params.m;
                Ciphertext rot_inner = automorph(inner, g, keys, ring);
                result = add(result, rot_inner, ring);
            }
        }
    }
    double bsgs_s2c = phase.sec();
    cout << " " << fixed << setprecision(2) << bsgs_s2c << " s" << endl;

    // BSGS C2S (same cost structure)
    cout << "  [2/3] DigitExtraction (PS, degree " << (params.p-1) << ")..."
         << " [extrapolated] " << fixed << setprecision(0) << de_time << " s" << endl;

    cout << "  [3/3] CoeffToSlot (BSGS)..." << flush;
    phase.start();
    {
        HoistedDecomp hd = hoist(ct, keys, ring);
        vector<Ciphertext> baby_rots(baby, Ciphertext(ct.modulus, ct.ptxt_space));
        baby_rots[0] = ct;
        for (long b = 1; b < baby; b++) {
            long g = 1;
            for (long i = 0; i < b; i++) g = (g * params.g_ell) % params.m;
            baby_rots[b] = apply_rotation(ct, hd, g, keys, ring);
        }
        Ciphertext result(ct.modulus, ct.ptxt_space);
        NTL::clear(result.c0); NTL::clear(result.c1);
        long w_inv = InvMod(65529L, params.p);
        long D_inv = InvMod((long)params.D, params.p);
        for (long gs = 0; gs < giant; gs++) {
            Ciphertext inner(ct.modulus, ct.ptxt_space);
            NTL::clear(inner.c0); NTL::clear(inner.c1);
            for (long b = 0; b < baby; b++) {
                long d = b + gs * baby;
                if (d >= params.D) break;
                vector<long> diag_vals(params.D);
                for (long i = 0; i < params.D; i++) {
                    long j = (i + d) % params.D;
                    long val = PowerMod(w_inv, (long)(((long long)i*j) % params.D), params.p);
                    diag_vals[i] = (long)(((long long)val * D_inv) % params.p);
                }
                ZZX diag_poly; ring.encode_diagonal(diag_poly, diag_vals, params.p);
                Ciphertext term = mul_const(baby_rots[b], diag_poly, ring);
                inner = add(inner, term, ring);
            }
            if (gs == 0) { result = inner; }
            else {
                long g = 1; long step = gs*baby;
                for (long i = 0; i < step; i++) g = (g * params.g_ell) % params.m;
                Ciphertext rot_inner = automorph(inner, g, keys, ring);
                result = add(result, rot_inner, ring);
            }
        }
    }
    double bsgs_c2s = phase.sec();
    cout << " " << fixed << setprecision(2) << bsgs_c2s << " s" << endl;

    double baseline_linear = bsgs_s2c + bsgs_c2s;
    double baseline_total = baseline_linear + de_time;
    cout << "\n  BASELINE TOTAL:" << endl;
    cout << "    Linear (S2C+C2S): " << fixed << setprecision(2) << baseline_linear << " s" << endl;
    cout << "    Digit Extraction: " << fixed << setprecision(0) << de_time << " s" << endl;
    cout << "    Bootstrap TOTAL:  " << baseline_total << " s ("
         << baseline_total/60 << " min)" << endl;

    // ---- TOWER: Optimized Bootstrap ----
    cout << "\n================================================================" << endl;
    cout << "  TOWER-OPTIMIZED Bootstrap (butterfly linear transform)" << endl;
    cout << "================================================================" << endl;

    TowerTransform tower(params, ring);

    cout << "\n  [1/3] SlotToCoeff (Tower, " << tower.layers.size() << " layers)..." << flush;
    phase.start();
    Ciphertext tower_s2c_ct = tower.apply(ct, keys);
    double tower_s2c = phase.sec();
    cout << " " << fixed << setprecision(2) << tower_s2c << " s" << endl;

    cout << "  [2/3] DigitExtraction (PS, same as baseline)..."
         << " [extrapolated] " << fixed << setprecision(0) << de_time << " s" << endl;

    cout << "  [3/3] CoeffToSlot (Tower)..." << flush;
    phase.start();
    Ciphertext tower_c2s_ct = tower.apply_inverse(ct, keys);
    double tower_c2s = phase.sec();
    cout << " " << fixed << setprecision(2) << tower_c2s << " s" << endl;

    double tower_linear = tower_s2c + tower_c2s;
    double tower_total = tower_linear + de_time;
    cout << "\n  TOWER TOTAL:" << endl;
    cout << "    Linear (S2C+C2S): " << fixed << setprecision(2) << tower_linear << " s" << endl;
    cout << "    Digit Extraction: " << fixed << setprecision(0) << de_time << " s" << endl;
    cout << "    Bootstrap TOTAL:  " << tower_total << " s ("
         << tower_total/60 << " min)" << endl;

    // ---- Final Comparison ----
    cout << "\n================================================================" << endl;
    cout << "  COMPARISON" << endl;
    cout << "================================================================" << endl;
    cout << fixed << setprecision(2);
    cout << "\n  +---------------------+-----------+-----------+---------+" << endl;
    cout << "  | Phase               | Baseline  | Tower     | Speedup |" << endl;
    cout << "  +---------------------+-----------+-----------+---------+" << endl;
    cout << "  | SlotToCoeff         | " << setw(7) << bsgs_s2c << " s | "
         << setw(7) << tower_s2c << " s | "
         << setw(5) << bsgs_s2c/tower_s2c << "x |" << endl;
    cout << "  | Digit Extraction    | " << setw(7) << (long)de_time << " s | "
         << setw(7) << (long)de_time << " s |  1.00x |" << endl;
    cout << "  | CoeffToSlot         | " << setw(7) << bsgs_c2s << " s | "
         << setw(7) << tower_c2s << " s | "
         << setw(5) << bsgs_c2s/tower_c2s << "x |" << endl;
    cout << "  +---------------------+-----------+-----------+---------+" << endl;
    cout << "  | Linear transform    | " << setw(7) << baseline_linear << " s | "
         << setw(7) << tower_linear << " s | "
         << setw(5) << baseline_linear/tower_linear << "x |" << endl;
    cout << "  | TOTAL Bootstrap     | " << setw(7) << (long)baseline_total << " s | "
         << setw(7) << (long)tower_total << " s | "
         << setw(5) << baseline_total/tower_total << "x |" << endl;
    cout << "  +---------------------+-----------+-----------+---------+" << endl;

    cout << "\n  Key findings:" << endl;
    cout << "    - Linear transform speedup: " << baseline_linear/tower_linear << "x" << endl;
    cout << "    - Total bootstrap speedup:  " << baseline_total/tower_total << "x" << endl;
    cout << "    - Digit extraction dominates: "
         << fixed << setprecision(1) << 100.0*de_time/baseline_total << "% of baseline, "
         << 100.0*de_time/tower_total << "% of tower" << endl;

    cout << "\n  Wall time: " << fixed << setprecision(1) << timer.sec() << " s" << endl;
    cout << "================================================================" << endl;
    return 0;
}
