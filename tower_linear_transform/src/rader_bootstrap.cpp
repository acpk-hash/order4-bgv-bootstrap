/*
 * rader_bootstrap.cpp - Rader-based linear transform verification
 * Computes all Rader stages and verifies product == Vandermonde
 * Uses HElib's internal G (slot polynomial) via template function
 */
#include <helib/helib.h>
#include <helib/EvalMap.h>
#include <helib/matmul.h>
#include <NTL/lzz_pX.h>
#include <NTL/lzz_pE.h>
#include <iostream>
#include <chrono>
#include <vector>

using namespace helib;
using namespace NTL;
using namespace std;

struct Timer {
    chrono::high_resolution_clock::time_point t0;
    void start() { t0 = chrono::high_resolution_clock::now(); }
    double sec() { return chrono::duration<double>(chrono::high_resolution_clock::now()-t0).count(); }
};

// Template function that has access to HElib's slot polynomial G
template <typename type>
static bool verifyRaderDecomposition(const EncryptedArray& _ea, long ell, long g_ell, long cofactor) {
    PA_INJECT(type)
    const long D = ell - 1;  // 96
    
    const EncryptedArrayDerived<type>& ea = _ea.getDerived(type());
    RBak bak; bak.save();
    _ea.getAlMod().restoreContext();
    const RX& G = ea.getG();
    long p_val = _ea.getAlMod().getZMStar().getP();
    long d = deg(G);
    
    cerr << "  G degree = " << d << endl;
    
    // Set up extension field with HElib's G
    REBak ebak; ebak.save();
    ea.restoreContextForG();
    
    // Compute omega_D (primitive D-th root of unity in F_{p^d})
    ZZ pd = power(ZZ(p_val), d);
    ZZ expD = (pd - 1) / D;
    
    RE omega_D;
    for (long trial = 2; trial < 200; trial++) {
        RX gen_poly;
        SetCoeff(gen_poly, 0, trial);
        SetCoeff(gen_poly, 1, 1);
        RE gen = conv<RE>(gen_poly);
        omega_D = power(gen, expD);
        // Verify order is exactly D
        if (!IsOne(omega_D) && IsOne(power(omega_D, D))) {
            // Check it's not a smaller divisor
            bool good = true;
            if (D % 2 == 0 && IsOne(power(omega_D, D/2))) good = false;
            if (D % 3 == 0 && IsOne(power(omega_D, D/3))) good = false;
            if (good) break;
        }
    }
    cerr << "  omega_" << D << " found (order " << D << " verified: " 
         << IsOne(power(omega_D, D)) << ")" << endl;
    
    // alpha = X^cofactor mod G (HElib's evaluation point base)
    RX X_poly;
    SetCoeff(X_poly, 1, 1);
    RE X_elem = conv<RE>(X_poly);
    RE alpha = power(X_elem, cofactor);
    cerr << "  alpha = X^" << cofactor << " mod G" << endl;
    
    // Build evaluation points: point[j] = alpha^{g^j mod ell}
    Vec<RE> points;
    points.SetLength(D);
    for (long j = 0; j < D; j++) {
        long exp = 1;
        for (long i = 0; i < j; i++) exp = (exp * g_ell) % ell;
        points[j] = power(alpha, exp);
    }
    
    // Build actual Vandermonde: V[i][j] = point[j]^i
    Mat<RX> V_actual;
    V_actual.SetDims(D, D);
    for (long j = 0; j < D; j++) NTL::set(V_actual[0][j]);
    for (long i = 1; i < D; i++)
        for (long j = 0; j < D; j++) {
            RE prod = conv<RE>(V_actual[i-1][j]) * points[j];
            V_actual[i][j] = rep(prod);
        }
    
    // Build Rader kernel: b[t] = alpha^{g^t mod ell} = point[t]
    // NTT(b): B[k] = sum_t b[t] * omega_D^{tk}
    Vec<RE> B;
    B.SetLength(D);
    for (long k = 0; k < D; k++) {
        clear(B[k]);
        for (long t = 0; t < D; t++) {
            long exp = (long)(((long long)t * k) % D);
            RE tw = power(omega_D, exp);
            B[k] += points[t] * tw;
        }
    }
    
    // Build Rader product: V_rader = P_out * INTT * diag(B) * conj_NTT
    // conj_NTT[k][j] = omega_D^{-jk}
    // INTT[s][k] = (1/D) * omega_D^{-sk}
    // diag(B)[k][k] = B[k]
    // P_out[i][dlog[i]] = 1 for i=1,...,D-1
    
    // Compute dlog table
    vector<long> dlog(ell, -1);
    long g_pow = 1;
    for (long s = 0; s < D; s++) {
        dlog[g_pow] = s;
        g_pow = (g_pow * g_ell) % ell;
    }
    
    // Compute product directly: V_rader[i][j] for i=1,...,D-1
    // V_rader[i][j] = (INTT * diag(B) * conj_NTT)[dlog[i]][j]
    // = (1/D) * sum_k omega_D^{-dlog[i]*k} * B[k] * omega_D^{-jk}
    // = (1/D) * sum_k B[k] * omega_D^{-(dlog[i]+j)*k}
    
    RE omega_inv = power(omega_D, D - 1);
    long D_inv_scalar = InvMod(D, p_val);
    RX D_inv_poly;
    SetCoeff(D_inv_poly, 0, D_inv_scalar);
    RE D_inv_E = conv<RE>(D_inv_poly);
    
    Mat<RX> V_rader;
    V_rader.SetDims(D, D);
    
    // Row 0: all ones
    for (long j = 0; j < D; j++) NTL::set(V_rader[0][j]);
    
    // Rows 1,...,D-1
    for (long i = 1; i < D; i++) {
        long s = dlog[i];  // i = g^s mod ell
        for (long j = 0; j < D; j++) {
            // V_rader[i][j] = (1/D) * sum_k B[k] * omega_D^{-(s+j)*k}
            RE val; clear(val);
            for (long k = 0; k < D; k++) {
                long exp = (long)((long long)(s + j) * k % D);
                RE tw = power(omega_inv, exp);
                val += B[k] * tw;
            }
            val *= D_inv_E;
            V_rader[i][j] = rep(val);
        }
    }
    
    // Compare
    bool match = true;
    long first_i = -1, first_j = -1;
    for (long i = 0; i < D && match; i++)
        for (long j = 0; j < D && match; j++)
            if (V_rader[i][j] != V_actual[i][j]) {
                match = false;
                first_i = i; first_j = j;
            }
    
    cerr << "  Rader product == Vandermonde? " << (match ? "YES ✓" : "NO ✗") << endl;
    if (!match) {
        cerr << "  First mismatch at [" << first_i << "][" << first_j << "]" << endl;
    }
    
    return match;
}

int main() {
    Timer timer;
    long p = 65537, m = 50731, r = 1, bits = 1500, c = 3;
    long ell = 97, g_ell = 5;
    vector<long> gens = {48117, 5239};
    vector<long> ords = {96, 29};
    vector<long> mvec_v = {97, 523};

    cout << "=== Rader Decomposition Verification ===" << endl;
    cout << "  m=" << m << ", p=" << p << ", ell=" << ell << ", D=" << (ell-1) << endl;

    cout << "[1] Context..." << flush;
    timer.start();
    Context context = ContextBuilder<BGV>()
        .m(m).p(p).r(r).bits(bits).c(c)
        .gens(gens).ords(ords).mvec(mvec_v)
        .buildModChain(true)
        .build();
    cout << " " << timer.sec() << "s" << endl;

    const EncryptedArray& ea = context.getEA();
    long cofactor = m / ell;  // 523

    cout << "[2] Verifying Rader decomposition..." << endl;
    timer.start();
    bool ok = verifyRaderDecomposition<PA_zz_p>(ea, ell, g_ell, cofactor);
    cout << "  Time: " << timer.sec() << "s" << endl;

    cout << "\n=== Result: " << (ok ? "PASS" : "FAIL") << " ===" << endl;
    return ok ? 0 : 1;
}
