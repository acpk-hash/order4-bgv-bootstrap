#pragma once
#include "params.h"
#include "poly_ring.h"
#include "ciphertext.h"
#include "keys.h"
#include <vector>

namespace tower_bgv {

// Tower Transform: evaluates the slot-to-coefficient linear transform
// using O(log D) layers, each with 2-3 diagonals.
//
// For D=96=32×3:
//   Layer 1 (stride 32): radix-3 butterfly (3 diags)
//   Layers 2-6 (strides 16,8,4,2,1): radix-2 butterfly (2 diags each)
//   Total: 6 hoistings + 12 rotations
//
// With layer merging (pairs of layers share hoisting):
//   3 merged groups: (stride 32+16), (stride 8+4), (stride 2+1)
//   Total: 3 hoistings + 12 rotations

struct TowerLayer {
    long stride;
    long degree;  // 2 or 3
    // Twiddle factors for this layer (one per butterfly group)
    // For degree-2: out[j] = in[j] + tw[j]*in[j+stride]
    //               out[j+s] = in[j] - tw[j]*in[j+stride]
    // These are PLAINTEXT constants (encoded as polynomials in the slot ring)
    std::vector<NTL::ZZX> twiddles;
    long automorph_k;  // the ring automorphism X->X^k for this stride
};

class TowerTransform {
public:
    const Params& params;
    const PolyRing& ring;
    std::vector<TowerLayer> layers;
    long omega_D;  // primitive D-th root of unity mod p (0 if doesn't exist)

    TowerTransform(const Params& p, const PolyRing& r);

    // Apply the full tower transform to a ciphertext (Slot-to-Coeff)
    Ciphertext apply(const Ciphertext& ct, const Keys& keys) const;

    // Apply inverse tower transform (Coeff-to-Slot)
    Ciphertext apply_inverse(const Ciphertext& ct, const Keys& keys) const;

    // Apply a single layer (for debugging)
    Ciphertext apply_layer(const Ciphertext& ct, long layer_idx,
                           const Keys& keys) const;

    // Count total operations
    long count_hoistings() const;
    long count_rotations() const;

private:
    void build_twiddles();
};

} // namespace tower_bgv
