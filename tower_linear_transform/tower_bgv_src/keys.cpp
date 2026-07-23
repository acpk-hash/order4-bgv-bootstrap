#include "keys.h"
#include <NTL/ZZX.h>
#include <NTL/ZZ.h>
#include <iostream>

namespace tower_bgv {

void Keys::gen_secret(const PolyRing& ring, long hamming_weight) {
    hwt = hamming_weight;
    ring.sparse_poly(secret, hwt);
}

void Keys::gen_ks_matrix(long k, const PolyRing& ring) {
    KSMatrix ks;
    ks.from_auto = k;
    ks.num_digits = 3;  // use 3 digits
    // digit_base B = 2^{ceil(log2(Q)/num_digits)} so that B^num_digits >= Q
    NTL::ZZ B;
    long bits = NTL::NumBits(modulus);
    long B_bits = (bits + ks.num_digits - 1) / ks.num_digits;
    B = NTL::ZZ(1) << B_bits;
    ks.digit_base = B;

    // s_k = s(X^k) mod Phi_m
    NTL::ZZX s_k;
    ring.automorph(s_k, secret, k);

    ks.b.resize(ks.num_digits);
    ks.a.resize(ks.num_digits);

    for (long i = 0; i < ks.num_digits; i++) {
        // a_i <- random
        ring.random_poly(ks.a[i], NTL::to_long(modulus/2));

        // e_i <- small error
        NTL::ZZX e_i;
        ring.random_poly(e_i, 1);

        // b_i = -a_i * s + p*e_i + B^i * s_k (mod Q, mod Phi_m)
        NTL::ZZX as_i;
        ring.mul_mod(as_i, ks.a[i], secret, modulus);

        NTL::ZZ Bi = NTL::power(B, i);
        NTL::ZZX Bi_sk;
        for (long j = 0; j <= NTL::deg(s_k); j++)
            NTL::SetCoeff(Bi_sk, j, NTL::coeff(s_k, j) * Bi);

        NTL::ZZX pe_i;
        for (long j = 0; j <= NTL::deg(e_i); j++)
            NTL::SetCoeff(pe_i, j, NTL::coeff(e_i, j) * p);

        ks.b[i] = Bi_sk + pe_i - as_i;
        balanced_reduce(ks.b[i], modulus);
    }

    ks_matrices[k] = std::move(ks);
}

void Keys::gen_relin_key(const PolyRing& ring) {
    // Relinearization key encrypts s^2 under s
    // Same structure as KS matrix: b_i + s*a_i ≈ B^i * s^2 (mod Q)
    relin_key.from_auto = 0;
    relin_key.num_digits = 3;
    long bits = NTL::NumBits(modulus);
    long B_bits = (bits + relin_key.num_digits - 1) / relin_key.num_digits;
    NTL::ZZ B = NTL::ZZ(1) << B_bits;
    relin_key.digit_base = B;

    // s^2 mod Phi_m
    NTL::ZZX s2;
    NTL::MulMod(s2, secret, secret, ring.phi_m);

    relin_key.b.resize(relin_key.num_digits);
    relin_key.a.resize(relin_key.num_digits);

    for (long i = 0; i < relin_key.num_digits; i++) {
        ring.random_poly(relin_key.a[i], NTL::to_long(modulus/2));

        NTL::ZZX e_i;
        ring.random_poly(e_i, 1);

        NTL::ZZX as_i;
        ring.mul_mod(as_i, relin_key.a[i], secret, modulus);

        NTL::ZZ Bi = NTL::power(B, i);
        NTL::ZZX Bi_s2;
        for (long j = 0; j <= NTL::deg(s2); j++)
            NTL::SetCoeff(Bi_s2, j, NTL::coeff(s2, j) * Bi);

        NTL::ZZX pe_i;
        for (long j = 0; j <= NTL::deg(e_i); j++)
            NTL::SetCoeff(pe_i, j, NTL::coeff(e_i, j) * p);

        relin_key.b[i] = Bi_s2 + pe_i - as_i;
        balanced_reduce(relin_key.b[i], modulus);
    }
    has_relin_key = true;
}

void Keys::gen_tower_ks(const Params& params, const PolyRing& ring) {
    std::cerr << "Generating tower KS matrices..." << std::endl;
    // For each tower stride, generate KS for σ^stride and σ^{-stride}
    for (long stride : params.tower_strides) {
        long k_fwd = stride;  // σ^stride in additive notation
        // In multiplicative notation: the automorphism is X -> X^{g^stride}
        // where g is the primitive root mod ell
        // For the Galois group action on slots: σ^k maps slot j to slot j+k
        // The actual ring automorphism is X -> X^{g^k mod m}
        long g_power = 1;
        for (long i = 0; i < stride; i++)
            g_power = (g_power * params.g_ell) % params.m;

        gen_ks_matrix(g_power, ring);
        std::cerr << "  KS for stride " << stride
                  << " (X->X^" << g_power << ") generated" << std::endl;

        // Inverse: σ^{-stride}
        long g_inv = 1;
        long g_ell_inv = NTL::InvMod(params.g_ell, params.m);
        for (long i = 0; i < stride; i++)
            g_inv = (g_inv * g_ell_inv) % params.m;

        gen_ks_matrix(g_inv, ring);
        std::cerr << "  KS for stride -" << stride
                  << " (X->X^" << g_inv << ") generated" << std::endl;
    }
    std::cerr << "  Total: " << ks_matrices.size() << " KS matrices" << std::endl;
}

// Digit decomposition: c = sum c_i * B^i where |c_i| < B/2
static std::vector<NTL::ZZX> digit_decompose(const NTL::ZZX& c, long num_digits,
                                              const NTL::ZZ& B) {
    std::vector<NTL::ZZX> digits(num_digits);
    NTL::ZZX remainder = c;
    NTL::ZZ half_B = B / 2;

    for (long d = 0; d < num_digits; d++) {
        for (long i = 0; i <= NTL::deg(remainder); i++) {
            NTL::ZZ coef = NTL::coeff(remainder, i);
            NTL::ZZ digit = coef % B;
            if (digit > half_B) digit -= B;
            if (digit < -half_B) digit += B;
            NTL::SetCoeff(digits[d], i, digit);
            NTL::SetCoeff(remainder, i, (coef - digit) / B);
        }
        remainder.normalize();
        digits[d].normalize();
    }
    return digits;
}

Ciphertext automorph(const Ciphertext& ct, long k,
                     const Keys& keys, const PolyRing& ring) {
    Ciphertext result(ct.modulus, ct.ptxt_space);

    // Step 1: Apply automorphism to both parts
    NTL::ZZX c0_k, c1_k;
    ring.automorph(c0_k, ct.c0, k);
    ring.automorph(c1_k, ct.c1, k);
    balanced_reduce(c0_k, ct.modulus);
    balanced_reduce(c1_k, ct.modulus);

    // Step 2: Key-switch c1_k from s(X^k) to s(X)
    auto it = keys.ks_matrices.find(k);
    if (it == keys.ks_matrices.end()) {
        std::cerr << "ERROR: No KS matrix for automorphism " << k << std::endl;
        result.c0 = c0_k;
        result.c1 = c1_k;
        return result;
    }
    const KSMatrix& ks = it->second;

    // Digit decompose c1_k
    std::vector<NTL::ZZX> digits = digit_decompose(c1_k, ks.num_digits, ks.digit_base);

    // result.c0 = c0_k + sum_i digits[i] * ks.b[i]
    // result.c1 = sum_i digits[i] * ks.a[i]
    result.c0 = c0_k;
    NTL::clear(result.c1);

    for (long i = 0; i < ks.num_digits; i++) {
        NTL::ZZX db, da;
        ring.mul_mod(db, digits[i], ks.b[i], ct.modulus);
        ring.mul_mod(da, digits[i], ks.a[i], ct.modulus);
        result.c0 += db;
        result.c1 += da;
    }
    balanced_reduce(result.c0, ct.modulus);
    balanced_reduce(result.c1, ct.modulus);

    result.noise_bound = ct.noise_bound * 2;  // rough estimate
    return result;
}

Ciphertext relinearize(const Ciphertext3& ct3, const Keys& keys, const PolyRing& ring) {
    Ciphertext result(ct3.modulus, ct3.ptxt_space);
    result.c0 = ct3.c0;
    result.c1 = ct3.c1;

    const KSMatrix& rk = keys.relin_key;
    std::vector<NTL::ZZX> digits = digit_decompose(ct3.c2, rk.num_digits, rk.digit_base);

    for (long i = 0; i < rk.num_digits; i++) {
        NTL::ZZX db, da;
        ring.mul_mod(db, digits[i], rk.b[i], ct3.modulus);
        ring.mul_mod(da, digits[i], rk.a[i], ct3.modulus);
        result.c0 += db;
        result.c1 += da;
    }
    balanced_reduce(result.c0, ct3.modulus);
    balanced_reduce(result.c1, ct3.modulus);

    result.noise_bound = ct3.noise_bound * 2;
    return result;
}

HoistedDecomp hoist(const Ciphertext& ct, const Keys& keys, const PolyRing& ring) {
    HoistedDecomp hd;
    // Use the first KS matrix's digit base for decomposition
    NTL::ZZ B = keys.ks_matrices.begin()->second.digit_base;
    long num_digits = keys.ks_matrices.begin()->second.num_digits;
    hd.digits = digit_decompose(ct.c1, num_digits, B);
    return hd;
}

Ciphertext apply_rotation(const Ciphertext& ct, const HoistedDecomp& hd,
                          long k, const Keys& keys, const PolyRing& ring) {
    Ciphertext result(ct.modulus, ct.ptxt_space);

    // Apply automorphism to c0
    ring.automorph(result.c0, ct.c0, k);
    balanced_reduce(result.c0, ct.modulus);

    // Key-switch using pre-decomposed digits
    auto it = keys.ks_matrices.find(k);
    if (it == keys.ks_matrices.end()) {
        std::cerr << "ERROR: No KS matrix for " << k << std::endl;
        return result;
    }
    const KSMatrix& ks = it->second;

    NTL::clear(result.c1);
    for (long i = 0; i < ks.num_digits; i++) {
        // Apply automorphism to each digit, then multiply by KS hint
        NTL::ZZX digit_k;
        ring.automorph(digit_k, hd.digits[i], k);
        balanced_reduce(digit_k, ct.modulus);

        NTL::ZZX db, da;
        ring.mul_mod(db, digit_k, ks.b[i], ct.modulus);
        ring.mul_mod(da, digit_k, ks.a[i], ct.modulus);
        result.c0 += db;
        result.c1 += da;
    }
    balanced_reduce(result.c0, ct.modulus);
    balanced_reduce(result.c1, ct.modulus);

    result.noise_bound = ct.noise_bound * 2;
    return result;
}

} // namespace tower_bgv
