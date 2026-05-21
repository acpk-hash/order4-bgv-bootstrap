#pragma once
#include "params.h"
#include "poly_ring.h"
#include <NTL/ZZX.h>
#include <NTL/ZZ.h>
#include <vector>

namespace tower_bgv {

// BGV Ciphertext: (c0, c1) such that c0 + s*c1 = m (mod p) (mod Q)
struct Ciphertext {
    NTL::ZZX c0, c1;
    NTL::ZZ modulus;      // current ciphertext modulus Q
    long ptxt_space;      // plaintext space p
    double noise_bound;   // estimated noise

    Ciphertext() : ptxt_space(0), noise_bound(0) {}
    Ciphertext(const NTL::ZZ& mod, long p) : modulus(mod), ptxt_space(p), noise_bound(0) {}
};

// Reduce all coefficients of f to [-Q/2, Q/2]
void balanced_reduce(NTL::ZZX& f, const NTL::ZZ& Q);

// Encrypt plaintext m under secret key s
Ciphertext encrypt(const NTL::ZZX& m, const NTL::ZZX& s,
                   const PolyRing& ring, const NTL::ZZ& Q, long p);

// Decrypt ciphertext under secret key s
NTL::ZZX decrypt(const Ciphertext& ct, const NTL::ZZX& s, const PolyRing& ring);

// Add two ciphertexts
Ciphertext add(const Ciphertext& a, const Ciphertext& b, const PolyRing& ring);

// Multiply ciphertext by plaintext constant
Ciphertext mul_const(const Ciphertext& ct, const NTL::ZZX& c, const PolyRing& ring);

// Degree-2 ciphertext from multiplication: c0 + s*c1 + s^2*c2 = m_a * m_b (mod p)
struct Ciphertext3 {
    NTL::ZZX c0, c1, c2;
    NTL::ZZ modulus;
    long ptxt_space;
    double noise_bound;

    Ciphertext3() : ptxt_space(0), noise_bound(0) {}
    Ciphertext3(const NTL::ZZ& mod, long p) : modulus(mod), ptxt_space(p), noise_bound(0) {}
};

// Multiply two ciphertexts (produces degree-2 ciphertext)
Ciphertext3 mul(const Ciphertext& a, const Ciphertext& b, const PolyRing& ring);

} // namespace tower_bgv
