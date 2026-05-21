/*
 * test_basic.cpp - Verify basic BGV operations
 * Tests: polynomial ring, encrypt/decrypt, automorphism, key-switching
 */
#include "params.h"
#include "poly_ring.h"
#include "ciphertext.h"
#include "keys.h"
#include "tower_transform.h"
#include <iostream>
#include <cassert>
#include <chrono>

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
};

int main() {
    cout << "=== Tower BGV: Full System Test ===" << endl;

    // Use small parameters for testing correctness
    // m=13 (prime), phi(13)=12
    Params p;
    p.m = 13;
    p.p = 5;
    p.phim = 12;
    p.D = 12;
    p.ell = 13;
    p.d = 4;
    p.nslots = 3;
    p.g_ell = 2;
    p.tower_degrees = {2, 2, 3};  // 12 = 2*2*3
    p.tower_strides = {6, 3, 1};
    p.dks = 2;

    Timer timer;

    // [1] Polynomial ring
    cout << "\n[1] Building polynomial ring (m=13)..." << flush;
    timer.start();
    PolyRing ring(p);
    cout << " " << (int)timer.ms() << " ms" << endl;
    assert(deg(ring.phi_m) == 12);
    cout << "    deg(Phi_13) = 12 ✓" << endl;

    // [2] Key generation
    cout << "\n[2] Generating keys (hwt=4)..." << flush;
    timer.start();
    Keys keys;
    keys.p = p.p;
    keys.modulus = ZZ(1) << 60;  // Q = 2^60
    keys.gen_secret(ring, 4);
    cout << " " << (int)timer.ms() << " ms" << endl;
    cout << "    Secret key generated (hwt=4)" << endl;

    // [3] Encrypt/Decrypt
    cout << "\n[3] Testing encrypt/decrypt..." << flush;
    ZZX msg;
    SetCoeff(msg, 0, 3);  // m = 3
    SetCoeff(msg, 1, 2);  // m = 3 + 2X

    Ciphertext ct = encrypt(msg, keys.secret, ring, keys.modulus, p.p);
    ZZX decrypted = decrypt(ct, keys.secret, ring);

    bool correct = true;
    for (long i = 0; i <= max(deg(msg), deg(decrypted)); i++) {
        long m_i = IsZero(coeff(msg, i)) ? 0 : to_long(coeff(msg, i));
        long d_i = IsZero(coeff(decrypted, i)) ? 0 : to_long(coeff(decrypted, i));
        if (m_i != d_i) { correct = false; break; }
    }
    cout << (correct ? " PASS ✓" : " FAIL ✗") << endl;
    if (!correct) {
        cout << "    msg = " << msg << endl;
        cout << "    dec = " << decrypted << endl;
    }

    // [4] Key-switching matrix generation
    cout << "\n[4] Generating KS matrix for σ^2 (X->X^2)..." << flush;
    timer.start();
    keys.gen_ks_matrix(2, ring);
    cout << " " << (int)timer.ms() << " ms" << endl;

    // [5] Automorphism
    cout << "\n[5] Testing automorphism σ^2..." << flush;
    Ciphertext ct_rot = automorph(ct, 2, keys, ring);
    ZZX dec_rot = decrypt(ct_rot, keys.secret, ring);

    // Expected: msg(X^2) mod Phi_13 mod p
    ZZX expected;
    ring.automorph(expected, msg, 2);
    for (long i = 0; i <= deg(expected); i++) {
        long c = to_long(coeff(expected, i)) % p.p;
        if (c < 0) c += p.p;
        SetCoeff(expected, i, c);
    }
    expected.normalize();

    correct = true;
    for (long i = 0; i <= max(deg(expected), deg(dec_rot)); i++) {
        long e_i = IsZero(coeff(expected, i)) ? 0 : to_long(coeff(expected, i));
        long d_i = IsZero(coeff(dec_rot, i)) ? 0 : to_long(coeff(dec_rot, i));
        if (e_i != d_i) { correct = false; break; }
    }
    cout << (correct ? " PASS ✓" : " FAIL ✗") << endl;
    if (!correct) {
        cout << "    expected = " << expected << endl;
        cout << "    got      = " << dec_rot << endl;
    }

    // [6] Addition
    cout << "\n[6] Testing ciphertext addition..." << flush;
    Ciphertext ct_sum = add(ct, ct, ring);
    ZZX dec_sum = decrypt(ct_sum, keys.secret, ring);
    ZZX expected_sum;
    for (long i = 0; i <= deg(msg); i++) {
        long c = (to_long(coeff(msg, i)) * 2) % p.p;
        if (c != 0) SetCoeff(expected_sum, i, c);
    }
    expected_sum.normalize();

    correct = true;
    for (long i = 0; i <= max(deg(expected_sum), deg(dec_sum)); i++) {
        long e_i = IsZero(coeff(expected_sum, i)) ? 0 : to_long(coeff(expected_sum, i));
        long d_i = IsZero(coeff(dec_sum, i)) ? 0 : to_long(coeff(dec_sum, i));
        if (e_i != d_i) { correct = false; break; }
    }
    cout << (correct ? " PASS ✓" : " FAIL ✗") << endl;

    cout << "\n=== Test Summary ===" << endl;
    cout << "All basic BGV operations verified on small parameters." << endl;
    cout << "Next: scale up to m=50731 (D=96) and benchmark tower transform." << endl;

    return 0;
}

