// Fails unless NTL is configured as on the machine that produced the paper's
// numbers: 60-bit single-precision primes, threads, GMP.
#include <NTL/ZZ.h>
#include <NTL/version.h>
#include <cstdio>
int main() {
  bool ok = NTL_SP_NBITS == 60;
#ifdef NTL_ENABLE_AVX_FFT
  ok = false;
#endif
#if !defined(NTL_THREADS) || !defined(NTL_GMP_LIP)
  ok = false;
#endif
  std::printf("NTL %s, NTL_SP_NBITS=%d: %s\n", NTL_VERSION, (int)NTL_SP_NBITS,
              ok ? "ok" : "MISMATCH with the measurement host");
  return ok ? 0 : 1;
}
