// combined_bench.cpp
// Measures the COMBINED operating point of [18]'s rank-matching rule with our
// order-four factored digit extraction at p = 257:
//   the factored polynomial Q (deg 55) only needs slot rank d >= 56, so the
//   [18]-compatible ring can shrink from M = 2^16 (d = 256, needed by their
//   deg-257 f and by unfactored P_A, deg 223) down to M = 2^15 (d = 128) or
//   M = 2^14 (d = 64). Same 128 slots in all cases; ring dimension N drops
//   32768 -> 16384 -> 8192. We evaluate the SAME factored path (y = x^4,
//   Q(y), 129*x + x^3*Q(y)) with the same plain BSGS evaluator on each ring.
// Usage: combined_bench <N>   with N in {8192, 16384, 32768}.
#include <iostream>
#include <vector>
#include <chrono>
#include <cmath>
#include <cassert>
#include "seal/seal.h"

using namespace seal;
using namespace std;

static size_t kswitch_counter = 0;
static inline void log_relin() { ++kswitch_counter; }

static const std::vector<uint64_t> Q_COEFFS = {91,65,245,23,253,38,251,198,170,166,23,194,31,88,233,48,233,152,153,7,31,153,34,157,83,206,154,60,207,99,106,230,196,218,246,155,100,76,86,146,81,74,206,184,163,48,228,16,112,9,53,175,252,72,23,224};
static const uint64_t C1 = 129;

static void bsgs_eval(const std::vector<uint64_t>& coeffs, const Ciphertext& x_enc,
                      const Evaluator& eval, const RelinKeys& rk,
                      const Modulus& plain_mod, Ciphertext& dst) {
    size_t deg = coeffs.size() - 1;
    while (deg > 0 && coeffs[deg] == 0) --deg;
    size_t bs = (size_t)std::ceil(std::sqrt((double)(deg + 1)));
    if (bs < 1) bs = 1;
    size_t gs = (deg / bs) + 1;
    auto mkplain = [&](uint64_t v) { Plaintext p; p.resize(1); p[0] = plain_mod.reduce(v); return p; };
    std::cerr << "    [bsgs] deg=" << deg << " bs=" << bs << " gs=" << gs << "\n";
    std::vector<Ciphertext> powers(bs);
    powers[0] = x_enc;
    for (size_t i = 1; i < bs; ++i) {
        eval.multiply(powers[i-1], x_enc, powers[i]);
        eval.relinearize_inplace(powers[i], rk); log_relin();
    }
    Ciphertext G = powers[bs-1];
    std::vector<Ciphertext> Gp(gs);
    if (gs > 1) Gp[1] = G;
    for (size_t k = 2; k < gs; ++k) {
        eval.multiply(Gp[k-1], G, Gp[k]);
        eval.relinearize_inplace(Gp[k], rk); log_relin();
    }
    Ciphertext acc; bool acc_set = false;
    for (size_t k = 0; k < gs; ++k) {
        Ciphertext block; bool block_set = false;
        uint64_t const_term = 0; bool has_const = false;
        for (size_t i = 0; i < bs; ++i) {
            size_t idx = k * bs + i;
            if (idx > deg) break;
            uint64_t c = coeffs[idx] % plain_mod.value();
            if (c == 0) continue;
            if (i == 0) { const_term = c; has_const = true; continue; }
            Ciphertext term;
            eval.multiply_plain(powers[i-1], mkplain(c), term);
            if (!block_set) { block = term; block_set = true; }
            else eval.add_inplace(block, term);
        }
        if (block_set && has_const) eval.add_plain_inplace(block, mkplain(const_term));
        if (!block_set && !has_const) continue;
        Ciphertext piece;
        if (k == 0) {
            if (block_set) piece = block;
            else { eval.multiply_plain(x_enc, mkplain(0), piece); eval.add_plain_inplace(piece, mkplain(const_term)); }
        } else {
            if (!block_set) { eval.multiply_plain(Gp[k], mkplain(const_term), piece); }
            else { eval.multiply(block, Gp[k], piece); eval.relinearize_inplace(piece, rk); log_relin(); }
        }
        if (!acc_set) { acc = piece; acc_set = true; }
        else eval.add_inplace(acc, piece);
    }
    dst = acc;
}

int main(int argc, char** argv) {
    size_t N = (argc > 1) ? (size_t)atol(argv[1]) : 32768;
    size_t d = (N == 8192) ? 64 : (N == 16384) ? 128 : 256;

    EncryptionParameters parms(scheme_type::bfv);
    parms.set_poly_modulus_degree(N);
    bool reduced_sec = false;
    if (N == 8192) {
        // Q path needs ~360 bits of budget; the 128-bit-security default at N=8192
        // (218 bits) is too small, so we inflate the modulus. Microbenchmark only,
        // reduced lattice security (disclosed), same convention as the D-trend rings.
        parms.set_coeff_modulus(CoeffModulus::Create(N, {50,50,50,50,50,50,50,44}));
        reduced_sec = true;
    } else {
        parms.set_coeff_modulus(CoeffModulus::BFVDefault(N));
    }
    parms.set_plain_modulus(257);
    SEALContext context(parms, true, reduced_sec ? sec_level_type::none : sec_level_type::tc128);
    const Modulus& plain_mod = context.first_context_data()->parms().plain_modulus();
    std::cout << "p=257 ring: N=" << N << " (M=2^" << (int)std::log2(2*N) << "), slot rank d=" << d
              << ", slots n=" << (N / d) << (reduced_sec ? "  [inflated modulus, reduced security]" : "  [128-bit security default]")
              << std::endl;
    std::cout << "rank check: Q needs d >= 56 -> " << (d >= 56 ? "OK" : "FAIL")
              << " ; their f (deg 257) needs d >= 257 -> " << (d >= 257 ? "OK" : "NOT AVAILABLE on this ring")
              << " ; unfactored P_A (deg 223) needs d >= 224 -> " << (d >= 224 ? "OK" : "NOT AVAILABLE on this ring")
              << std::endl;

    KeyGenerator keygen(context);
    SecretKey sk = keygen.secret_key();
    PublicKey pk; keygen.create_public_key(pk);
    RelinKeys rk; keygen.create_relin_keys(rk);
    Encryptor encryptor(context, pk);
    Evaluator evaluator(context);
    Decryptor decryptor(context, sk);

    auto encrypt_const = [&](long c) {
        Plaintext xp; xp.resize(1); xp[0] = (uint64_t)((c % 257 + 257) % 257);
        Ciphertext ce; encryptor.encrypt(xp, ce); return ce;
    };
    auto const_result = [&](const Ciphertext& ct) {
        Plaintext pt; decryptor.decrypt(ct, pt);
        return (long)(pt.coeff_count() ? pt[0] : 0);
    };

    long x_bench = (5*16 + 6) % 257;   // eta=5, lambda=6 -> x=86
    Ciphertext x_enc = encrypt_const(x_bench);
    std::cout << "encrypted constant input x=" << x_bench << "; noise budget "
              << decryptor.invariant_noise_budget(x_enc) << " bits" << std::endl;

    // probes
    { Ciphertext sq; evaluator.square(x_enc, sq); evaluator.relinearize_inplace(sq, rk);
      std::cout << "PROBE square(86): " << const_result(sq) << " expect " << (86*86 % 257) << std::endl; }

    // factored path timing at this ring
    auto run_factored = [&](const Ciphertext& xin, size_t& ks_out, long& budget_out) {
        kswitch_counter = 0;
        Ciphertext x2, x3, x4;
        evaluator.square(xin, x2);       evaluator.relinearize_inplace(x2, rk); log_relin();
        evaluator.multiply(xin, x2, x3); evaluator.relinearize_inplace(x3, rk); log_relin();
        evaluator.square(x2, x4);        evaluator.relinearize_inplace(x4, rk); log_relin();
        Ciphertext qy; bsgs_eval(Q_COEFFS, x4, evaluator, rk, plain_mod, qy);
        Ciphertext term1; { Plaintext c1p; c1p.resize(1); c1p[0] = plain_mod.reduce(C1); evaluator.multiply_plain(xin, c1p, term1); }
        Ciphertext term2; evaluator.multiply(x3, qy, term2); evaluator.relinearize_inplace(term2, rk); log_relin();
        Ciphertext outc; evaluator.add(term1, term2, outc);
        ks_out = kswitch_counter;
        budget_out = decryptor.invariant_noise_budget(outc);
        return outc;
    };

    size_t ks; long budget;
    auto t0 = std::chrono::steady_clock::now();
    Ciphertext fac = run_factored(x_enc, ks, budget);
    auto t1 = std::chrono::steady_clock::now();
    long ms = std::chrono::duration_cast<std::chrono::milliseconds>(t1 - t0).count();
    std::cout << "[combined] factored Q path on N=" << N << ": " << ms << " ms, "
              << ks << " key-switches, noise budget " << budget << std::endl;
    std::cout << "result for x=86 (eta=5,lambda=6): " << const_result(fac) << " expect 6" << std::endl;

    // correctness on a spread of support points
    int okc = 0, tot = 0;
    for (auto [e, l] : std::vector<std::pair<int,int>>{{-7,-7},{-5,3},{0,0},{3,-2},{7,7}}) {
        long x = (((e*16 + l) % 257) + 257) % 257;
        Ciphertext xe = encrypt_const(x);
        size_t k2; long b2;
        Ciphertext r = run_factored(xe, k2, b2);
        long got = const_result(r), expl = ((l % 257) + 257) % 257;
        ++tot; if (got == expl) ++okc;
        std::cout << "  support (eta=" << e << ",lam=" << l << ",x=" << x << "): got " << got
                  << " expect " << expl << (got == expl ? "  OK" : "  MISMATCH") << std::endl;
    }
    std::cout << "correctness: " << okc << "/" << tot << " support points OK" << std::endl;
    return 0;
}
