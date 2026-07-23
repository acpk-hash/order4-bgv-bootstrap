#include "poly_eval.h"
#include <cmath>
#include <iostream>

namespace tower_bgv {

Ciphertext poly_eval(const Ciphertext& ct, const std::vector<NTL::ZZ>& coeffs,
                     const Keys& keys, const PolyRing& ring) {
    long d = coeffs.size() - 1;
    while (d > 0 && NTL::IsZero(coeffs[d])) d--;

    if (d == 0) {
        Ciphertext result(ct.modulus, ct.ptxt_space);
        NTL::SetCoeff(result.c0, 0, coeffs[0]);
        NTL::clear(result.c1);
        result.noise_bound = 0;
        return result;
    }

    // Paterson-Stockmeyer: split f(x) = sum_{i=0}^{m-1} f_i(x) * x^{i*k}
    // where k = ceil(sqrt(d+1)), m = ceil((d+1)/k)
    // f_i is a polynomial of degree < k
    long k = (long)std::ceil(std::sqrt((double)(d + 1)));
    long m = (d + k) / k;

    std::cerr << "  poly_eval: degree=" << d << ", k=" << k << ", m=" << m << std::endl;

    // Step 1: Precompute powers ct^1, ct^2, ..., ct^k
    std::vector<Ciphertext> powers(k + 1);
    powers[1] = ct;
    for (long i = 2; i <= k; i++) {
        Ciphertext3 prod = mul(powers[i-1], ct, ring);
        powers[i] = relinearize(prod, keys, ring);
        if (i % 50 == 0)
            std::cerr << "    power " << i << "/" << k << std::endl;
    }

    // Step 2: Evaluate using Horner on the "giant steps"
    // result = f_{m-1}(ct) * ct^{(m-1)*k} + ... + f_1(ct) * ct^k + f_0(ct)
    // Using Horner: result = (...((f_{m-1} * ct^k + f_{m-2}) * ct^k + f_{m-3}) ... ) * ct^k + f_0

    // Evaluate each f_i(ct) = sum_{j=0}^{k-1} coeffs[i*k+j] * ct^j
    auto eval_baby = [&](long block_idx) -> Ciphertext {
        Ciphertext result(ct.modulus, ct.ptxt_space);
        NTL::clear(result.c0);
        NTL::clear(result.c1);
        result.noise_bound = 0;

        long base = block_idx * k;
        // Start with constant term
        if (base <= d && !NTL::IsZero(coeffs[base])) {
            NTL::SetCoeff(result.c0, 0, coeffs[base]);
        }

        for (long j = 1; j < k && base + j <= d; j++) {
            if (NTL::IsZero(coeffs[base + j])) continue;
            // Add coeffs[base+j] * ct^j
            NTL::ZZX scalar;
            NTL::SetCoeff(scalar, 0, coeffs[base + j]);
            Ciphertext term = mul_const(powers[j], scalar, ring);
            result = add(result, term, ring);
        }
        return result;
    };

    // Horner's method on giant steps
    Ciphertext result = eval_baby(m - 1);
    for (long i = m - 2; i >= 0; i--) {
        // result = result * ct^k + f_i(ct)
        Ciphertext3 prod = mul(result, powers[k], ring);
        result = relinearize(prod, keys, ring);
        Ciphertext fi = eval_baby(i);
        result = add(result, fi, ring);

        if ((m - 2 - i) % 20 == 0)
            std::cerr << "    giant step " << (m - 2 - i) << "/" << (m - 1) << std::endl;
    }

    return result;
}

} // namespace tower_bgv
