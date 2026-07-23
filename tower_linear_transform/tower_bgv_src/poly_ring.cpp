#include "poly_ring.h"
#include <NTL/ZZXFactoring.h>
#include <NTL/BasicThreadPool.h>
#include <algorithm>
#include <random>

namespace tower_bgv {

// Compute cyclotomic polynomial Phi_m(X)
static NTL::ZZX compute_cyclotomic(long m) {
    // Use NTL's built-in cyclotomic polynomial computation
    // Phi_m(X) = prod_{d|m} (X^d - 1)^{mu(m/d)}
    NTL::ZZX phi;
    NTL::SetCoeff(phi, 0, 1); // start with 1

    // Simple method: Phi_m = (X^m - 1) / prod_{d|m, d<m} Phi_d
    // For prime m: Phi_m = 1 + X + X^2 + ... + X^{m-1}
    // For m = p*q (coprime): use factorization

    // Use the formula: Phi_m(X) = GCD of X^m-1 and appropriate factors
    // Actually, NTL can factor cyclotomic polynomials directly
    // Let's just build it from the definition for our specific m values

    if (m == 50731) {
        // m = 97 * 523
        // Phi_{50731}(X) has degree phi(50731) = phi(97)*phi(523) = 96*522 = 50112
        // This is too large to store explicitly!
        // Instead, we work with the SLOT polynomial (degree d=6 or 18)
        // and the dimension D=96

        // Actually for BGV, we work in Z[X]/Phi_m(X) where deg(Phi_m) = phi(m)
        // phi(50731) = 96 * 522 = 50112
        // This is the ring dimension N = 50112

        // Build Phi_m using the formula for m = p*q (coprime):
        // Phi_{pq}(X) = Phi_p(X^q) / Phi_p(X) ... no that's wrong
        // Phi_{pq}(X) = prod of (X^{pq/d} - 1)^{mu(d)} for d | pq

        // Simpler: for m=p*q with p,q prime:
        // Phi_{pq}(X) = (X^{pq}-1)(X-1) / ((X^p-1)(X^q-1))
        // = (1+X^q+X^{2q}+...+X^{(p-1)q}) evaluated at X -> sum

        // Actually the cleanest way:
        // Phi_{97*523}(X) = Phi_{97}(X^{523}) / gcd(Phi_{97}(X^{523}), Phi_{523}(X^{97}))
        // For coprime p,q: Phi_{pq}(X) = Phi_p(X^q) / Phi_p(X) when gcd(p,q)=1... no

        // Correct formula for squarefree m = p1*p2*...*pk:
        // Use inclusion-exclusion with Mobius function
        // For m = 97*523: Phi_m(X) = (X^m - 1)(X - 1) / ((X^97 - 1)(X^523 - 1))

        NTL::ZZX num, den, xm, x97, x523, x1;
        NTL::SetCoeff(xm, m, 1); NTL::SetCoeff(xm, 0, -1);   // X^m - 1
        NTL::SetCoeff(x97, 97, 1); NTL::SetCoeff(x97, 0, -1); // X^97 - 1
        NTL::SetCoeff(x523, 523, 1); NTL::SetCoeff(x523, 0, -1); // X^523 - 1
        NTL::SetCoeff(x1, 1, 1); NTL::SetCoeff(x1, 0, -1);   // X - 1

        // Phi_m = (X^m - 1)(X - 1) / ((X^97 - 1)(X^523 - 1))
        NTL::mul(num, xm, x1);
        NTL::mul(den, x97, x523);
        NTL::div(phi, num, den);
    } else {
        // Generic: compute using NTL
        // For small m, just use the definition
        NTL::ZZX xm_minus_1;
        NTL::SetCoeff(xm_minus_1, m, 1);
        NTL::SetCoeff(xm_minus_1, 0, -1);

        // Divide out all Phi_d for proper divisors d of m
        phi = xm_minus_1;
        for (long d = 1; d < m; d++) {
            if (m % d == 0) {
                NTL::ZZX phi_d = compute_cyclotomic(d);
                NTL::ZZX q, r;
                NTL::DivRem(q, r, phi, phi_d);
                if (NTL::IsZero(r))
                    phi = q;
            }
        }
    }

    return phi;
}

PolyRing::PolyRing(const Params& p) : params(p) {
    phi_m = compute_cyclotomic(p.m);
}

void PolyRing::reduce(NTL::ZZX& f) const {
    NTL::ZZX q, r;
    NTL::DivRem(q, r, f, phi_m);
    f = r;
}

void PolyRing::automorph(NTL::ZZX& result, const NTL::ZZX& f, long k) const {
    // f(X) -> f(X^k) mod Phi_m(X)
    long n = NTL::deg(f);
    NTL::clear(result);
    for (long i = 0; i <= n; i++) {
        if (!NTL::IsZero(NTL::coeff(f, i))) {
            long new_exp = (i * k) % params.m;
            // X^{new_exp} mod Phi_m(X) - for cyclotomic, X^m = 1 effectively
            // but we need proper reduction
            NTL::ZZX term;
            NTL::SetCoeff(term, new_exp, NTL::coeff(f, i));
            result += term;
        }
    }
    reduce(result);
}

void PolyRing::mul_mod(NTL::ZZX& result, const NTL::ZZX& a, const NTL::ZZX& b,
                       const NTL::ZZ& modulus) const {
    NTL::MulMod(result, a, b, phi_m);
    // Reduce coefficients mod modulus
    for (long i = 0; i <= NTL::deg(result); i++) {
        NTL::ZZ c = NTL::coeff(result, i) % modulus;
        if (c > modulus/2) c -= modulus;
        NTL::SetCoeff(result, i, c);
    }
    result.normalize();
}

void PolyRing::add_mod(NTL::ZZX& result, const NTL::ZZX& a, const NTL::ZZX& b,
                       const NTL::ZZ& modulus) const {
    result = a + b;
    for (long i = 0; i <= NTL::deg(result); i++) {
        NTL::ZZ c = NTL::coeff(result, i) % modulus;
        if (c > modulus/2) c -= modulus;
        NTL::SetCoeff(result, i, c);
    }
    result.normalize();
}

void PolyRing::random_poly(NTL::ZZX& f, long bound) const {
    NTL::clear(f);
    for (long i = 0; i < params.phim; i++) {
        long c = NTL::RandomBnd(2*bound + 1) - bound;
        if (c != 0) NTL::SetCoeff(f, i, c);
    }
}

void PolyRing::sparse_poly(NTL::ZZX& f, long h) const {
    NTL::clear(f);
    std::vector<long> positions(params.phim);
    for (long i = 0; i < params.phim; i++) positions[i] = i;

    std::mt19937 rng(std::random_device{}());
    std::shuffle(positions.begin(), positions.end(), rng);

    for (long i = 0; i < h; i++) {
        long sign = (NTL::RandomBnd(2) == 0) ? 1 : -1;
        NTL::SetCoeff(f, positions[i], sign);
    }
}

void PolyRing::compute_slot_roots(std::vector<long>& roots, long p) const {
    // Find primitive m-th root of unity mod p
    // Need: p ≡ 1 (mod m) so that ζ_m exists in Z_p
    long m = params.m;
    long g_p = 0;
    for (long g = 2; g < p; g++) {
        bool is_prim = true;
        // Check g is primitive root mod p by testing g^{(p-1)/q} != 1 for prime q | p-1
        long pm1 = p - 1;
        for (long q = 2; q * q <= pm1; q++) {
            if (pm1 % q == 0) {
                if (NTL::PowerMod(g, pm1 / q, p) == 1) { is_prim = false; break; }
                while (pm1 % q == 0) pm1 /= q;
            }
        }
        if (pm1 > 1 && NTL::PowerMod(g, (p-1) / pm1, p) == 1) is_prim = false;
        if (is_prim) { g_p = g; break; }
    }

    long zeta_m = NTL::PowerMod(g_p, (p - 1) / m, p);

    // Roots of Phi_m mod p, ordered by g_ell^j
    // root[j] = ζ_m^{g_ell^j mod m}
    roots.resize(params.D);
    for (long j = 0; j < params.D; j++) {
        long exp = 1;
        for (long i = 0; i < j; i++) exp = (exp * params.g_ell) % m;
        roots[j] = NTL::PowerMod(zeta_m, exp, p);
    }
}

void PolyRing::crt_encode(NTL::ZZX& out, const std::vector<long>& slots, long p) const {
    std::vector<long> roots;
    compute_slot_roots(roots, p);
    long D = params.D;

    // Lagrange interpolation: find f(X) of degree < D such that f(roots[j]) = slots[j]
    NTL::clear(out);
    for (long i = 0; i < D; i++) {
        if (slots[i] == 0) continue;
        // Compute Lagrange basis L_i(X) = prod_{j!=i} (X - roots[j]) / (roots[i] - roots[j])
        long denom = 1;
        for (long j = 0; j < D; j++) {
            if (j == i) continue;
            denom = (long)(((long long)denom * ((roots[i] - roots[j] + p) % p)) % p);
        }
        long inv_denom = NTL::InvMod(denom, p);

        // Build numerator polynomial prod_{j!=i} (X - roots[j])
        std::vector<long> num(D, 0);
        num[0] = 1;
        long deg_num = 0;
        for (long j = 0; j < D; j++) {
            if (j == i) continue;
            long neg_rj = (p - roots[j]) % p;
            for (long k = deg_num + 1; k >= 1; k--)
                num[k] = (long)(((long long)num[k-1] + (long long)num[k] * neg_rj) % p);
            num[0] = (long)((long long)num[0] * neg_rj % p);
            deg_num++;
        }

        // Add slots[i] * inv_denom * num to out
        long scale = (long)(((long long)slots[i] * inv_denom) % p);
        for (long k = 0; k < D; k++) {
            long c = (long)(((long long)scale * num[k]) % p);
            if (c != 0) {
                NTL::ZZ cur = NTL::coeff(out, k);
                NTL::SetCoeff(out, k, cur + NTL::ZZ(c));
            }
        }
    }
    // Reduce coefficients mod p
    for (long k = 0; k <= NTL::deg(out); k++) {
        long c = NTL::IsZero(NTL::coeff(out, k)) ? 0 : NTL::to_long(NTL::coeff(out, k)) % p;
        if (c < 0) c += p;
        NTL::SetCoeff(out, k, c);
    }
    out.normalize();
}

void PolyRing::crt_decode(std::vector<long>& slots, const NTL::ZZX& poly, long p) const {
    std::vector<long> roots;
    compute_slot_roots(roots, p);
    long D = params.D;
    slots.resize(D);

    for (long j = 0; j < D; j++) {
        // Evaluate poly at roots[j] mod p
        long val = 0;
        long power = 1;
        for (long i = 0; i <= NTL::deg(poly); i++) {
            long c = NTL::IsZero(NTL::coeff(poly, i)) ? 0 : NTL::to_long(NTL::coeff(poly, i)) % p;
            if (c < 0) c += p;
            val = (long)((val + (long long)c * power) % p);
            power = (long)(((long long)power * roots[j]) % p);
        }
        slots[j] = val;
    }
}

void PolyRing::encode_diagonal(NTL::ZZX& out, const std::vector<long>& vals, long p) const {
    crt_encode(out, vals, p);
}

Params Params::create_D96() {
    Params p;
    p.m = 50731;
    p.p = 65537;
    p.ell = 97;
    p.phim = 50112;  // phi(50731) = 96*522
    p.D = 96;
    p.d = 6;         // ord_97(65537) = 6 (actually 18 for full ring, 6 for this dim)
    p.nslots = 16;
    p.g_ell = 5;     // primitive root mod 97

    // Modulus chain: use ~60-bit primes, enough for linear transform
    // For proof-of-concept, use fewer primes
    p.num_primes = 10;
    p.dks = 3;       // 3 digits for key-switching

    // Tower structure for D=96 = 32*3
    p.tower_degrees = {3, 2, 2, 2, 2, 2};  // top-down: first radix-3, then radix-2
    p.tower_strides = {32, 16, 8, 4, 2, 1};

    return p;
}

Params Params::create_D256() {
    Params p;
    p.m = 49601;
    p.p = 65537;
    p.ell = 257;
    p.phim = 49152;  // phi(49601) = 256*192
    p.D = 256;
    p.d = 16;
    p.nslots = 16;
    p.g_ell = 3;     // primitive root mod 257

    p.num_primes = 12;
    p.dks = 3;

    // Tower structure for D=256 = 2^8
    p.tower_degrees = {2, 2, 2, 2, 2, 2, 2, 2};
    p.tower_strides = {128, 64, 32, 16, 8, 4, 2, 1};

    return p;
}

} // namespace tower_bgv
