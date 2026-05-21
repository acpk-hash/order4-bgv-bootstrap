#include "ciphertext.h"
#include <NTL/ZZX.h>
#include <NTL/ZZ.h>

namespace tower_bgv {

void balanced_reduce(NTL::ZZX& f, const NTL::ZZ& Q) {
    NTL::ZZ half_Q = Q / 2;
    for (long i = 0; i <= NTL::deg(f); i++) {
        NTL::ZZ c = NTL::coeff(f, i) % Q;
        if (c > half_Q) c -= Q;
        if (c < -half_Q) c += Q;
        NTL::SetCoeff(f, i, c);
    }
    f.normalize();
}

Ciphertext encrypt(const NTL::ZZX& m, const NTL::ZZX& s,
                   const PolyRing& ring, const NTL::ZZ& Q, long p) {
    Ciphertext ct(Q, p);

    // a <- random in R_Q
    NTL::ZZX a;
    ring.random_poly(a, NTL::to_long(Q/2));

    // e <- small error
    NTL::ZZX e;
    ring.random_poly(e, 1);  // ternary error

    // c1 = a
    ct.c1 = a;

    // c0 = -a*s + p*e + m (mod Q) (mod Phi_m)
    NTL::ZZX as;
    ring.mul_mod(as, a, s, Q);

    ct.c0 = m;
    NTL::ZZX pe;
    for (long i = 0; i <= NTL::deg(e); i++)
        NTL::SetCoeff(pe, i, NTL::coeff(e, i) * p);

    ct.c0 += pe;
    ct.c0 -= as;
    balanced_reduce(ct.c0, Q);

    ct.noise_bound = p * 2.0;  // rough estimate
    return ct;
}

NTL::ZZX decrypt(const Ciphertext& ct, const NTL::ZZX& s, const PolyRing& ring) {
    // m = (c0 + s*c1) mod Q mod p
    NTL::ZZX sc1;
    ring.mul_mod(sc1, ct.c1, s, ct.modulus);

    NTL::ZZX result = ct.c0 + sc1;
    balanced_reduce(result, ct.modulus);

    // Reduce mod p
    for (long i = 0; i <= NTL::deg(result); i++) {
        long c = NTL::to_long(NTL::coeff(result, i) % ct.ptxt_space);
        if (c < 0) c += ct.ptxt_space;
        NTL::SetCoeff(result, i, c);
    }
    result.normalize();
    return result;
}

Ciphertext add(const Ciphertext& a, const Ciphertext& b, const PolyRing& ring) {
    Ciphertext result(a.modulus, a.ptxt_space);
    result.c0 = a.c0 + b.c0;
    result.c1 = a.c1 + b.c1;
    balanced_reduce(result.c0, a.modulus);
    balanced_reduce(result.c1, a.modulus);
    result.noise_bound = a.noise_bound + b.noise_bound;
    return result;
}

Ciphertext mul_const(const Ciphertext& ct, const NTL::ZZX& c, const PolyRing& ring) {
    Ciphertext result(ct.modulus, ct.ptxt_space);
    ring.mul_mod(result.c0, ct.c0, c, ct.modulus);
    ring.mul_mod(result.c1, ct.c1, c, ct.modulus);
    result.noise_bound = ct.noise_bound * 10;  // rough
    return result;
}

Ciphertext3 mul(const Ciphertext& a, const Ciphertext& b, const PolyRing& ring) {
    Ciphertext3 result(a.modulus, a.ptxt_space);
    // c0 = a.c0 * b.c0
    ring.mul_mod(result.c0, a.c0, b.c0, a.modulus);
    // c1 = a.c0 * b.c1 + a.c1 * b.c0
    NTL::ZZX t1, t2;
    ring.mul_mod(t1, a.c0, b.c1, a.modulus);
    ring.mul_mod(t2, a.c1, b.c0, a.modulus);
    ring.add_mod(result.c1, t1, t2, a.modulus);
    // c2 = a.c1 * b.c1
    ring.mul_mod(result.c2, a.c1, b.c1, a.modulus);
    result.noise_bound = a.noise_bound * b.noise_bound;
    return result;
}

} // namespace tower_bgv
