#pragma once
#include "params.h"
#include "poly_ring.h"
#include "ciphertext.h"
#include <NTL/ZZX.h>
#include <NTL/ZZ.h>
#include <vector>
#include <unordered_map>

namespace tower_bgv {

// Key-switch matrix for automorphism σ^k
// Encrypts s(X^k) under s(X): KS = {(b_i, a_i)} where b_i + s*a_i ≈ B^i * s(X^k)
struct KSMatrix {
    long from_auto;  // the automorphism index k
    std::vector<NTL::ZZX> b;  // KS hint top row
    std::vector<NTL::ZZX> a;  // KS hint bottom row (or use seed)
    long num_digits;
    NTL::ZZ digit_base;  // B
};

// Secret key + key-switching matrices
struct Keys {
    NTL::ZZX secret;       // secret key polynomial s
    long hwt;              // Hamming weight of s
    NTL::ZZ modulus;       // ciphertext modulus Q
    long p;                // plaintext modulus

    // Key-switching matrices indexed by automorphism power
    std::unordered_map<long, KSMatrix> ks_matrices;

    // Relinearization key: encrypts s^2 under s
    KSMatrix relin_key;
    bool has_relin_key = false;

    // Generate secret key with given Hamming weight
    void gen_secret(const PolyRing& ring, long hamming_weight);

    // Generate KS matrix for automorphism σ^k
    void gen_ks_matrix(long k, const PolyRing& ring);

    // Generate relinearization key (for s^2 -> s)
    void gen_relin_key(const PolyRing& ring);

    // Generate all KS matrices needed for tower decomposition
    void gen_tower_ks(const Params& params, const PolyRing& ring);
};

// Apply automorphism σ^k to ciphertext (includes key-switching)
Ciphertext automorph(const Ciphertext& ct, long k,
                     const Keys& keys, const PolyRing& ring);

// Relinearize: convert degree-2 ciphertext to degree-1
Ciphertext relinearize(const Ciphertext3& ct3, const Keys& keys, const PolyRing& ring);

// Hoisted automorphism: decompose once, apply multiple rotations
struct HoistedDecomp {
    std::vector<NTL::ZZX> digits;  // digit decomposition of c1
};

HoistedDecomp hoist(const Ciphertext& ct, const Keys& keys, const PolyRing& ring);

// Apply rotation using pre-hoisted decomposition
Ciphertext apply_rotation(const Ciphertext& ct, const HoistedDecomp& hd,
                          long k, const Keys& keys, const PolyRing& ring);

} // namespace tower_bgv
