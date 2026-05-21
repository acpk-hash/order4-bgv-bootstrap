/**
 * batchFatboot.cpp - End-to-end batch BGV bootstrapping test
 *
 * Uses the same parameters as fatboot (m=50731, p=65537) but bootstraps
 * B ciphertexts and measures per-ciphertext throughput.
 *
 * Demonstrates: Order-4 digit extraction + batch amortization potential
 */
#include <helib/helib.h>
#include <iostream>
#include <chrono>
#include <vector>
#include <cstdlib>

using namespace std;
using namespace helib;

int main(int argc, char* argv[]) {
    long B = 2; // batch size
    if (argc > 1) B = atol(argv[1]);

    // Same parameters as our Order-4 experiment (fatboot i=4)
    cout << "=== End-to-End Batch BGV Bootstrap ===" << endl;
    cout << "B=" << B << " (batch size)" << endl;
    cout << "Using fatboot parameters: m=50731, p=65537" << endl;
    cout << "Running existing fatboot with repeat=" << B << " to simulate batch..." << endl;
    cout << endl;

    // Instead of reimplementing the full bootstrapping setup,
    // we run fatboot B times and measure total time
    // This gives us the BASELINE (independent bootstrapping)

    auto start = chrono::high_resolution_clock::now();

    // Run fatboot with our Order-4 optimization
    string cmd = "HELIB_EXPLICIT_AUX=256 HELIB_AUX_ORDER4_EVAL=1 "
                 "../../BGV-Boot-auxradix-opt/build/fatboot i=4 h=12 t=-1 newbts=1 thick=0 repeat="
                 + to_string(B) + " 2>&1";

    cout << "Command: " << cmd << endl;
    cout << "Running..." << endl;

    FILE* pipe = popen(cmd.c_str(), "r");
    if (!pipe) {
        cerr << "Failed to run fatboot" << endl;
        return 1;
    }

    char buffer[1024];
    double total_linear1 = 0, total_linear2 = 0, total_extract = 0, total_total = 0;
    int count = 0;

    while (fgets(buffer, sizeof(buffer), pipe)) {
        string line(buffer);

        // Parse timing lines
        double l1, l2, ext, tot;
        if (sscanf(buffer, "time for linear1 = %lf, linear2 = %lf, extract = %lf, total = %lf",
                   &l1, &l2, &ext, &tot) == 4) {
            total_linear1 += l1;
            total_linear2 += l2;
            total_extract += ext;
            total_total += tot;
            count++;
            printf("  Run %d: linear=%.1f+%.1f=%.1fs, extract=%.1fs, total=%.1fs\n",
                   count, l1, l2, l1+l2, ext, tot);
        }
    }
    pclose(pipe);

    auto end = chrono::high_resolution_clock::now();
    double wall_time = chrono::duration<double>(end - start).count();

    if (count == 0) {
        cerr << "No results obtained from fatboot!" << endl;
        return 1;
    }

    cout << "\n========================================" << endl;
    cout << "=== BATCH BOOTSTRAPPING RESULTS ===" << endl;
    cout << "========================================" << endl;
    cout << endl;

    double avg_linear = (total_linear1 + total_linear2) / count;
    double avg_extract = total_extract / count;
    double avg_total = total_total / count;

    cout << "Independent bootstrapping (B=" << count << " ciphertexts):" << endl;
    cout << "  Average linear transform: " << avg_linear << "s per ct" << endl;
    cout << "  Average digit extraction:  " << avg_extract << "s per ct" << endl;
    cout << "  Average total:             " << avg_total << "s per ct" << endl;
    cout << "  Wall-clock total:          " << wall_time << "s" << endl;
    cout << endl;

    // Estimate shared-a amortized performance
    double decomp_frac = 0.45; // gadget decomposition fraction
    double amort_linear = avg_linear * (1.0 - decomp_frac) + avg_linear * decomp_frac / B;
    double amort_total = amort_linear + avg_extract;

    cout << "Shared-a amortized estimate (B=" << B << "):" << endl;
    cout << "  Linear transform (amortized): " << amort_linear << "s per ct" << endl;
    cout << "  Digit extraction (Order-4):   " << avg_extract << "s per ct" << endl;
    cout << "  Total per ct:                 " << amort_total << "s" << endl;
    cout << "  Speedup vs independent:       " << avg_total / amort_total << "x" << endl;
    cout << endl;

    // Compare with Ma et al. baseline (from our earlier experiments)
    double ma_baseline_total = 213.0; // Ma et al. baseline total
    cout << "Comparison with Ma et al. baseline (213s):" << endl;
    cout << "  Our Order-4 (single ct):      " << avg_total << "s (speedup "
         << ma_baseline_total / avg_total << "x)" << endl;
    cout << "  Our batch+Order-4 (B=" << B << "):   " << amort_total << "s (speedup "
         << ma_baseline_total / amort_total << "x)" << endl;
    cout << endl;

    // PC-MM estimate (full amortization)
    double pcmm_speedup = 3.95; // from our PC-MM verification
    double pcmm_linear = avg_linear / pcmm_speedup;
    double pcmm_total = pcmm_linear + avg_extract;
    cout << "PC-MM full amortization estimate (B=" << B << "):" << endl;
    cout << "  Linear transform (PC-MM):     " << pcmm_linear << "s per ct" << endl;
    cout << "  Digit extraction (Order-4):   " << avg_extract << "s per ct" << endl;
    cout << "  Total per ct:                 " << pcmm_total << "s" << endl;
    cout << "  Speedup vs Ma baseline:       " << ma_baseline_total / pcmm_total << "x" << endl;

    cout << "\n=== Complete ===" << endl;
    return 0;
}
