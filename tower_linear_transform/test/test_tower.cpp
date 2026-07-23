/*
 * test_tower.cpp - Tower transform correctness verification
 * Uses m=13, p=157 where D=12 | (p-1)=156, so NTT works over Z_p.
 */
#include "params.h"
#include "poly_ring.h"
#include "ciphertext.h"
#include "keys.h"
#include "tower_transform.h"
#include <iostream>
#include <cassert>
#include <vector>

using namespace tower_bgv;
using namespace NTL;
using namespace std;

int main() {
    cout << "=== Tower Transform Correctness Test ===" << endl;
    cout << "Parameters: m=13, p=157, D=12, d=1, nslots=12" << endl;
    cout << "Tower: radices [2,2,3], strides [6,3,1]" << endl;

    Params p;
    p.m = 13;
    p.p = 157;
    p.phim = 12;
    p.D = 12;
    p.ell = 13;
    p.d = 1;
    p.nslots = 12;
    p.g_ell = 2;
    p.tower_degrees = {2, 2, 3};
    p.tower_strides = {6, 3, 1};
    p.dks = 3;

    int pass = 0, fail = 0;

    // Setup
    cout << "\n[Setup] Building ring and keys..." << flush;
    PolyRing ring(p);
    assert(deg(ring.phi_m) == 12);

    Keys keys;
    keys.p = p.p;
    keys.modulus = ZZ(1) << 200;  // large enough for multiple tower stages
    keys.gen_secret(ring, 4);
    keys.gen_relin_key(ring);

    // Generate KS matrices for ALL rotations (needed for full DFT)
    long g_ell_inv = InvMod(p.g_ell, p.m);
    for (long stride = 1; stride < p.D; stride++) {
        long g_fwd = 1;
        for (long i = 0; i < stride; i++) g_fwd = (g_fwd * p.g_ell) % p.m;
        keys.gen_ks_matrix(g_fwd, ring);
    }
    cout << " done (" << keys.ks_matrices.size() << " KS matrices)" << endl;

    // [1] Test CRT encode/decode
    cout << "\n[1] Testing CRT encode/decode..." << flush;
    vector<long> slots = {3, 7, 2, 11, 5, 9, 1, 4, 8, 6, 10, 12};
    ZZX encoded;
    ring.crt_encode(encoded, slots, p.p);
    vector<long> decoded;
    ring.crt_decode(decoded, encoded, p.p);

    bool crt_ok = (slots == decoded);
    cout << (crt_ok ? " PASS" : " FAIL") << endl;
    if (!crt_ok) {
        cout << "    slots:   "; for (auto v : slots) cout << v << " "; cout << endl;
        cout << "    decoded: "; for (auto v : decoded) cout << v << " "; cout << endl;
    }
    if (crt_ok) pass++; else fail++;

    // [2] Test encrypt/decrypt preserves slots
    cout << "\n[2] Testing encrypt/decrypt with CRT..." << flush;
    Ciphertext ct = encrypt(encoded, keys.secret, ring, keys.modulus, p.p);
    ZZX dec_poly = decrypt(ct, keys.secret, ring);
    vector<long> dec_slots;
    ring.crt_decode(dec_slots, dec_poly, p.p);

    bool enc_ok = (slots == dec_slots);
    cout << (enc_ok ? " PASS" : " FAIL") << endl;
    if (!enc_ok) {
        cout << "    expected: "; for (auto v : slots) cout << v << " "; cout << endl;
        cout << "    got:      "; for (auto v : dec_slots) cout << v << " "; cout << endl;
    }
    if (enc_ok) pass++; else fail++;

    // [3] Test forward tower transform (butterfly NTT)
    cout << "\n[3] Testing forward tower transform..." << flush;
    TowerTransform tower(p, ring);

    Ciphertext ct_ntt = tower.apply(ct, keys);
    ZZX ntt_poly = decrypt(ct_ntt, keys.secret, ring);
    vector<long> ntt_slots;
    ring.crt_decode(ntt_slots, ntt_poly, p.p);

    // The butterfly NTT computes M3*M2*M1*x (DFT with digit-reversed input)
    // Verify it's non-trivial (not all zeros or same as input)
    bool nontrivial = (ntt_slots != slots);
    bool nonzero = false;
    for (long v : ntt_slots) if (v != 0) { nonzero = true; break; }
    cout << (nontrivial && nonzero ? " PASS (non-trivial output)" : " FAIL") << endl;
    if (nontrivial && nonzero) pass++; else fail++;
    cout << "    output: "; for (auto v : ntt_slots) cout << v << " "; cout << endl;

    // [4] Test round-trip: inverse(forward(x)) = x
    cout << "\n[4] Testing round-trip (inverse ∘ forward = identity)..." << flush;
    Ciphertext ct_inv = tower.apply_inverse(ct_ntt, keys);
    ZZX inv_poly = decrypt(ct_inv, keys.secret, ring);
    vector<long> inv_slots;
    ring.crt_decode(inv_slots, inv_poly, p.p);

    bool inv_ok = (slots == inv_slots);
    cout << (inv_ok ? " PASS" : " FAIL") << endl;
    if (!inv_ok) {
        cout << "    original: "; for (auto v : slots) cout << v << " "; cout << endl;
        cout << "    got:      "; for (auto v : inv_slots) cout << v << " "; cout << endl;
    }
    if (inv_ok) pass++; else fail++;

    // Summary
    cout << "\n=== Summary: " << pass << " passed, " << fail << " failed ===" << endl;
    return fail > 0 ? 1 : 0;
}
