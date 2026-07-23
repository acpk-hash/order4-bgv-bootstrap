#pragma once
#include "ciphertext.h"
#include "keys.h"
#include "poly_ring.h"
#include <NTL/ZZX.h>
#include <NTL/ZZ.h>
#include <vector>

namespace tower_bgv {

// Homomorphic polynomial evaluation using Paterson-Stockmeyer method
// Evaluates f(ct) where f is a plaintext polynomial of degree d
// Cost: O(sqrt(d)) multiplications + O(sqrt(d)) precomputed powers
Ciphertext poly_eval(const Ciphertext& ct, const std::vector<NTL::ZZ>& coeffs,
                     const Keys& keys, const PolyRing& ring);

} // namespace tower_bgv
