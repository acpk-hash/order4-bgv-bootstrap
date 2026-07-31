// parameter.h
//
// Recommended BGV bootstrapping parameter sets, revised against the attack of
// ePrint 2026/279 and against the stock Lattice Estimator.
//
// Security convention. Every set carries two LWE instances and its security is
// the minimum of the two:
//   main key         n = phi(m), the full modulus chain Q, weight h_main = 120,
//                    which is HElib's default and the value the prior work fixes
//   encapsulated key same n, a short bootstrapping modulus q_boot, weight h_enc
// The symbols follow ePrint 2024/115: h' would be the main key there. We write
// h_main and h_enc to avoid the clash with ePrint 2024/164, which swaps them.
//
// Scope of the attack. ePrint 2026/279 needs the ring to be Z_q[X]/(X^N + 1)
// with N a power of two. Phi_m(X) = X^{phi(m)} + 1 exactly when m is a power of
// two, so the test is whether the cyclotomic index is a power of two, and it
// admits no borderline case. General cyclotomic sets are therefore reported
// under the stock estimator alone; power-of-two sets carry both columns.
//
// lambda_stock  minimum over the stock estimator's six attacks
// lambda_after  minimum over every attack run including the rotated hybrid;
//               this is NOT lambda_stock minus a fixed delta, because the stock
//               minimum is often already cheaper than the rotated tool's own
//               plain hybrid
//
// log2Q is the chain HElib actually builds, not the value passed as `bits`.
// The two differ by roughly 300 bits at m = 2^16 and 210 at m = 2^17, and the
// offset is not a constant of the ring alone, so `bits_request` below is the
// value that produced the stated `log2Q` on our build and should be verified
// against the printed chain rather than assumed.

#pragma once
#include <cstddef>

namespace order4 {

struct ParamSet {
  const char* name;
  long        p;             // plaintext prime
  long        m;             // cyclotomic index
  long        phim;          // ring dimension, n
  long        A;             // order-four radix modulo p
  long        h_main;        // main key Hamming weight
  long        h_enc;         // encapsulated bootstrapping key weight
  long        bits_request;  // value to pass to HElib
  double      log2Q;         // chain it built on our machine
  long        slots;
  double      lambda_stock;  // stock estimator minimum
  double      lambda_after;  // minimum over all attacks, rotated included
  bool        power_of_two;  // true when the attack of 2026/279 is in scope
};

// ---------------------------------------------------------------------------
// Power-of-two rings. In scope of ePrint 2026/279.
//
// The shipped configuration sat at 1314 to 1332 bits and measures 73.6 to 74.8
// after the attack. An exhaustive modulus search at weight 120, in steps of 5
// bits, puts the 80 bit crossing at log2Q = 1200 (80.04, lost at 1205). We do
// not operate at the boundary: the sets below sit well inside it, because the
// attack model is recent and its heuristics may still evolve.
//
// The main-key weight is deliberately left at 120. Raising it does not work:
// once the chain is re-sized for the heavier key, as the sqrt(h) noise law
// requires, weight 270 reaches only 79.07 and the curve is flat from 330. The
// optimal guessing dimension never approaches zero, so the attack is never
// switched off, only made gradually less profitable.
// ---------------------------------------------------------------------------
static const ParamSet PO2_SETS[] = {
  // 80 bit target. Measured end to end for all ten plaintext primes below.
  {"po2-80",  65537, 65536,  32768, 256, 120, 26,  782, 1066.0, 32768, 95.6, 86.7, true},
  // 100 bit target. Larger ring; no shipped set exists at this ring.
  {"po2-100", 65537, 131072, 65536, 256, 120, 26, 1468, 1676.0, 32768, 117.8, 106.0, true},
  // 128 bit target. Reachable only by doubling the ring: at m = 2^16 the
  // modulus a 128 bit target permits is close to what bootstrapping alone
  // consumes, which leaves no usable capacity.
  {"po2-128", 65537, 131072, 65536, 256, 120, 26, 1069, 1277.0, 32768, 143.6, 131.4, true},
};

// The ten plaintext primes measured on m = 2^16 at the 80 bit set, with the
// order-four radix each one uses. A is not interchangeable between primes: it
// must satisfy A^2 = -1 mod p, and passing another prime's radix silently
// disables the order-four path.
struct Po2Prime { long p; long A; long slots; };
static const Po2Prime PO2_PRIMES[] = {
  { 1297,  36,     8}, { 1601,  40,    32}, { 2521,  71,     4},
  { 3137,  56,    32}, { 4513,  95,    16}, { 7057,  84,     8},
  {13457, 116,     8}, {14401, 120,    32}, {15377, 124,     8},
  {65537, 256, 32768},
};

// ---------------------------------------------------------------------------
// General cyclotomic rings. Out of scope of ePrint 2026/279, so lambda_after
// equals lambda_stock and only one column is meaningful.
//
// Chosen the way the prior work chooses: fix p, take m = q1*q2 with a small
// d = ord_m(p), size Q so the main key sits just above the target, then verify
// both instances. Our additional filter is p = 1 mod 4 together with the
// injectivity condition max(a,b) > 2B where p = a^2 + b^2, which caps h_enc per
// prime.
// ---------------------------------------------------------------------------
static const ParamSet GEN_SETS[] = {
  // 80 bit tier. Rings the paper already runs, with the chain lengthened until
  // the main key sits just above the target; the spare margin becomes capacity.
  {"gen-80-A",  4513, 37637, 36960,  95, 120, 12, 1101, 1510.0,  9240, 81.1, 81.1, false},
  {"gen-80-B", 14401, 43033, 41140, 120, 120, 12, 1228, 1680.0,  3740, 81.3, 81.3, false},
  {"gen-80-C", 13457, 45193, 44100, 116, 120, 12, 1318, 1800.0,  6300, 81.4, 81.4, false},
  {"gen-80-D",  2521, 50731, 50112,  71, 120, 12, 1498, 2040.0,  2784, 81.9, 81.9, false},
  {"gen-80-E", 65537, 50731, 50112, 256, 120, 12, 1500, 2017.0,  2784, 82.6, 82.6, false},
  // 128 bit tier. New instantiations of the same primes on larger rings.
  {"gen-128-A",13457, 65047, 62776, 116, 120, 26, 1039, 1427.0, 31388,128.5,128.5, false},
  {"gen-128-B", 2917, 53983, 52488,  54, 120, 26,  856, 1184.0, 13122,128.6,128.6, false},
  {"gen-128-C",15377, 55459, 53640, 124, 120, 26,  877, 1212.0, 13410,128.4,128.4, false},
  {"gen-128-D", 8101, 57521, 56832,  90, 120, 26,  934, 1287.0,  7104,128.5,128.5, false},
  {"gen-128-E",14401, 54359, 51480, 120, 120, 26,  839, 1161.0,  5148,128.5,128.5, false},
};

// ---------------------------------------------------------------------------
// Sets withdrawn as recommendations.
//
// These were taken from HElib's bundled parameter tables without a per-ring
// security check and sit at 56.9 to 76.1 bits under a fixed 1055 bit chain.
// Their timing measurements remain valid as demonstrations; they must not be
// used as parameters. Two of them, p = 2917 and p = 15377, fail on the
// encapsulated instance as well, so shortening the chain does not repair them.
// ---------------------------------------------------------------------------
static const ParamSet WITHDRAWN[] = {
  {"withdrawn-2917",   2917, 17587, 17136,  54, 120, 12, 0, 1055.0, 0, 56.9, 56.9, false},
  {"withdrawn-15377", 15377, 19679, 17880, 124, 120, 12, 0, 1055.0, 0, 59.4, 59.4, false},
  {"withdrawn-8101",   8101, 23771, 21600,  90, 120, 12, 0, 1055.0, 0, 69.5, 69.5, false},
  {"withdrawn-1601",   1601, 24901, 24192,  40, 120, 12, 0, 1055.0, 0, 76.1, 76.1, false},
};

static const size_t N_PO2       = sizeof(PO2_SETS) / sizeof(ParamSet);
static const size_t N_PO2_PRIME = sizeof(PO2_PRIMES) / sizeof(Po2Prime);
static const size_t N_GEN       = sizeof(GEN_SETS) / sizeof(ParamSet);
static const size_t N_WITHDRAWN = sizeof(WITHDRAWN) / sizeof(ParamSet);

}  // namespace order4
