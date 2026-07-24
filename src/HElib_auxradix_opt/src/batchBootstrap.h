#ifndef HELIB_BATCH_BOOTSTRAP_H
#define HELIB_BATCH_BOOTSTRAP_H

/**
 * @file batchBootstrap.h
 * @brief Batch BGV Bootstrapping with shared-a amortization
 *
 * Core idea: Multiple ciphertexts sharing the same 'a' component
 * can amortize the gadget decomposition cost during key-switching.
 * Combined with Order-4 character-filtered digit extraction for
 * accelerating the per-ciphertext nonlinear step.
 */

#include <helib/Ctxt.h>
#include <helib/Context.h>
#include <vector>

namespace helib {

/**
 * Batch automorphism on shared-a ciphertexts.
 *
 * Given B ciphertexts that share the same 'a' part (first component),
 * apply automorphism sigma^k to all of them, amortizing the gadget
 * decomposition of sigma^k(a) across all B ciphertexts.
 *
 * Standard cost: B * (decompose + inner_product) = B * O(N*l)
 * Amortized cost: 1 * decompose + B * inner_product = O(N*l) + B*O(N)
 * Savings: (B-1) * decompose_cost
 *
 * @param cts Vector of B ciphertexts (modified in place)
 * @param k The automorphism index (apply sigma^k)
 */
void batchSmartAutomorph(std::vector<Ctxt>& cts, long k);

/**
 * Batch linear transform (EvalMap) on shared-a ciphertexts.
 *
 * Applies the same linear transform to B ciphertexts, amortizing
 * key-switch costs across the batch.
 *
 * @param cts Vector of B ciphertexts (modified in place)
 * @param ea The EncryptedArray defining the linear transform
 * @param direction true for forward (SlotToCoeff), false for inverse
 */
void batchEvalMap(std::vector<Ctxt>& cts, const EncryptedArray& ea,
                  bool direction);

/**
 * Complete batch BGV bootstrapping.
 *
 * Combines:
 * 1. Batch linear transforms (shared-a amortization)
 * 2. Per-ciphertext Order-4 digit extraction
 *
 * @param cts Vector of B ciphertexts to bootstrap (modified in place)
 */
void batchThinReCrypt(std::vector<Ctxt>& cts);

} // namespace helib

#endif // HELIB_BATCH_BOOTSTRAP_H
