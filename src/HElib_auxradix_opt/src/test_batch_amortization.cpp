/**
 * test_batch_amortization.cpp
 * Measures key-switch amortization potential for batch BGV bootstrapping
 */
#include <helib/helib.h>
#include <iostream>
#include <chrono>
#include <vector>

using namespace std;
using namespace helib;

int main()
{
    long m = 4095;
    long p = 127;
    long r = 1;
    long bits = 300;
    long c = 2;

    cout << "=== Batch BGV Key-Switch Amortization Test ===" << endl;
    cout << "m=" << m << " p=" << p << " bits=" << bits << endl;

    Context context = ContextBuilder<BGV>()
        .m(m).p(p).r(r).bits(bits).c(c).build();

    long nslots = context.getNSlots();
    cout << "nslots=" << nslots << " security=" << context.securityLevel() << endl;

    SecKey secretKey(context);
    secretKey.GenSecKey();
    addSome1DMatrices(secretKey);
    const PubKey& pk = secretKey;
    const EncryptedArray& ea = context.getEA();

    // Create test ciphertexts
    vector<long> ptxt(nslots, 1);

    // Measure single automorphism cost
    long gen = context.getZMStar().ZmStarGen(0);
    cout << "Using automorphism generator: " << gen << endl;

    // Warm up
    {
        Ctxt tmp(pk);
        ea.encrypt(tmp, pk, ptxt);
        tmp.smartAutomorph(gen);
    }

    // Measure single key-switch time
    int trials = 5;
    double single_ks_time = 0;
    for (int t = 0; t < trials; t++) {
        Ctxt ct(pk);
        ea.encrypt(ct, pk, ptxt);
        auto start = chrono::high_resolution_clock::now();
        ct.smartAutomorph(gen);
        auto end = chrono::high_resolution_clock::now();
        single_ks_time += chrono::duration<double>(end - start).count();
    }
    single_ks_time /= trials;
    cout << "\nSingle automorphism time: " << single_ks_time*1000 << " ms" << endl;

    // Measure batch (B independent) automorphism time
    for (int B : {2, 4, 8}) {
        vector<Ctxt> cts;
        for (int i = 0; i < B; i++) {
            Ctxt ct(pk);
            ea.encrypt(ct, pk, ptxt);
            cts.push_back(ct);
        }

        double batch_time = 0;
        for (int t = 0; t < trials; t++) {
            vector<Ctxt> tmp_cts = cts;
            auto start = chrono::high_resolution_clock::now();
            for (int i = 0; i < B; i++) {
                tmp_cts[i].smartAutomorph(gen);
            }
            auto end = chrono::high_resolution_clock::now();
            batch_time += chrono::duration<double>(end - start).count();
        }
        batch_time /= trials;

        // Estimate shared-a savings
        // From HElib profiling: keySwitchPart has two phases:
        //   breakIntoDigits: ~40-50% of key-switch time (gadget decomposition)
        //   keySwitchDigits: ~50-60% (inner product with eval key)
        // With shared-a: breakIntoDigits done once, keySwitchDigits done B times
        double decomp_frac = 0.45; // from profiling
        double amortized_per_ct = single_ks_time * (1.0 - decomp_frac)
                                + single_ks_time * decomp_frac / B;
        double amortized_total = amortized_per_ct * B;

        cout << "\nB=" << B << ":" << endl;
        cout << "  Standard batch: " << batch_time*1000 << " ms total, "
             << batch_time/B*1000 << " ms/ct" << endl;
        cout << "  Shared-a (est): " << amortized_total*1000 << " ms total, "
             << amortized_per_ct*1000 << " ms/ct" << endl;
        cout << "  Speedup: " << batch_time / amortized_total << "x" << endl;
    }

    // Full bootstrapping estimate
    cout << "\n=== Full Bootstrapping Estimate ===" << endl;
    // From our profiling: linear transform has ~27 automorphisms
    // Each automorphism costs ~0.54s (from m=50731 experiment)
    // Digit extraction: 81s (with Order-4)
    double aut_cost = single_ks_time; // per automorphism
    long num_auts = 27;
    double T_lin = num_auts * aut_cost;
    double T_ext = 81.0; // Order-4 digit extraction (from experiment)

    cout << "Linear transform (single ct): " << T_lin << "s" << endl;
    cout << "Digit extraction (Order-4): " << T_ext << "s" << endl;
    cout << "Total (single ct): " << T_lin + T_ext << "s" << endl;

    for (int B : {2, 4, 8, 16}) {
        double decomp_frac = 0.45;
        double T_lin_batch = num_auts * (aut_cost * (1.0 - decomp_frac) * B
                                        + aut_cost * decomp_frac);
        double T_ext_batch = T_ext * B; // digit extraction not amortized
        double total_batch = T_lin_batch + T_ext_batch;
        double per_ct = total_batch / B;
        double baseline_per_ct = T_lin + T_ext;

        cout << "  B=" << B << ": per-ct=" << per_ct << "s"
             << " speedup=" << baseline_per_ct/per_ct << "x" << endl;
    }

    cout << "\n=== Done ===" << endl;
    return 0;
}
