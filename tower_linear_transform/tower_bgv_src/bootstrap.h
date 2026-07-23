#pragma once
#include "ciphertext.h"
#include "keys.h"
#include "poly_ring.h"
#include "tower_transform.h"
#include <NTL/ZZX.h>
#include <NTL/ZZ.h>
#include <vector>

namespace tower_bgv {

struct BootstrapTimings {
    double mod_raise_ms;
    double coeff_to_slot_ms;
    double digit_extract_ms;
    double slot_to_coeff_ms;
    double total_ms;
};

// Compute the digit extraction polynomial for BGV
// f(x) = x mod p, evaluated over Z_q
// Returns coefficients of the polynomial of degree <= q-1
std::vector<NTL::ZZ> digit_extract_poly(long p, long q);

// Full BGV bootstrap using tower transform
// Input: ciphertext at low modulus q_low
// Output: ciphertext at high modulus Q with refreshed noise
Ciphertext bootstrap(const Ciphertext& ct_low, const NTL::ZZ& Q_high,
                     const Keys& keys, const TowerTransform& tower,
                     const PolyRing& ring, BootstrapTimings& timings);

// Mod-raise: embed ciphertext from modulus q to modulus Q
Ciphertext mod_raise(const Ciphertext& ct, const NTL::ZZ& Q_new);

} // namespace tower_bgv
