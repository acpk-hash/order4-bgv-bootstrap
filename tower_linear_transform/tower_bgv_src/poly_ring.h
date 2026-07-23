#pragma once
#include "params.h"
#include <NTL/ZZX.h>
#include <NTL/ZZ_pX.h>
#include <NTL/ZZ_pXFactoring.h>
#include <NTL/ZZ.h>
#include <vector>

namespace tower_bgv {

// Polynomial ring R = Z[X] / Phi_m(X)
// All operations are mod Phi_m(X) and mod the current modulus Q
class PolyRing {
public:
    const Params& params;
    NTL::ZZX phi_m;           // cyclotomic polynomial Phi_m(X)
    NTL::ZZ_pXModulus phi_mod; // for fast modular reduction

    explicit PolyRing(const Params& p);

    // Reduce poly mod Phi_m(X)
    void reduce(NTL::ZZX& f) const;

    // Automorphism: f(X) -> f(X^k) mod Phi_m(X)
    void automorph(NTL::ZZX& result, const NTL::ZZX& f, long k) const;

    // Multiply mod Phi_m and mod current modulus
    void mul_mod(NTL::ZZX& result, const NTL::ZZX& a, const NTL::ZZX& b,
                 const NTL::ZZ& modulus) const;

    // Add mod modulus
    void add_mod(NTL::ZZX& result, const NTL::ZZX& a, const NTL::ZZX& b,
                 const NTL::ZZ& modulus) const;

    // Random polynomial with coefficients in [-bound, bound]
    void random_poly(NTL::ZZX& f, long bound) const;

    // Sparse random polynomial with exactly h non-zero ±1 coefficients
    void sparse_poly(NTL::ZZX& f, long h) const;

    // CRT slot operations (for d=1 case: slots are scalars in Z_p)
    // Compute roots of Phi_m mod p
    void compute_slot_roots(std::vector<long>& roots, long p) const;

    // Encode slot values into polynomial via Lagrange interpolation
    void crt_encode(NTL::ZZX& out, const std::vector<long>& slots, long p) const;

    // Decode polynomial to slot values (evaluate at roots)
    void crt_decode(std::vector<long>& slots, const NTL::ZZX& poly, long p) const;

    // Encode a diagonal: polynomial that puts val[j] in slot j
    void encode_diagonal(NTL::ZZX& out, const std::vector<long>& vals, long p) const;
};

} // namespace tower_bgv
