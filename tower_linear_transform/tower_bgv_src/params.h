#pragma once
#include <NTL/ZZX.h>
#include <NTL/ZZ_pX.h>
#include <NTL/lzz_pX.h>
#include <NTL/ZZ.h>
#include <vector>
#include <cstdint>

namespace tower_bgv {

struct Params {
    long m;          // cyclotomic index
    long p;          // plaintext modulus
    long phim;       // phi(m)
    long D;          // dimension of the linear transform (phi(ell))
    long ell;        // prime factor of m for the main dimension
    long d;          // slot polynomial degree = ord_ell(p)
    long nslots;     // number of slots = D/d

    // Modulus chain: Q = prod of primes q_i
    std::vector<long> primes;  // the primes in the chain
    long num_primes;           // number of primes (levels)

    // Key-switching parameters
    long dks;        // number of digits for key-switching

    // Tower structure
    std::vector<long> tower_degrees;  // [2,2,2,2,2,3] for D=96
    std::vector<long> tower_strides;  // [1,2,4,8,16,32] for D=96

    // Primitive roots
    long g_ell;      // primitive root mod ell (generator of (Z/ell)*)

    static Params create_D96();   // m=50731, D=96
    static Params create_D256();  // m=49601, D=256
};

} // namespace tower_bgv
