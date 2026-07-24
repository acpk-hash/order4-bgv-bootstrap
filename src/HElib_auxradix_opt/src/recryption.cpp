/* Copyright (C) 2012-2020 IBM Corp.
 * This program is Licensed under the Apache License, Version 2.0
 * (the "License"); you may not use this file except in compliance
 * with the License. You may obtain a copy of the License at
 *   http://www.apache.org/licenses/LICENSE-2.0
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License. See accompanying LICENSE file.
 */
#include <NTL/BasicThreadPool.h>

#include <helib/recryption.h>
#include <helib/EncryptedArray.h>
#include <helib/EvalMap.h>
#include <helib/powerful.h>
#include <helib/CtPtrs.h>
#include <helib/intraSlot.h>
#include <helib/norms.h>
#include <helib/sample.h>
#include <helib/debugging.h>
#include <helib/fhe_stats.h>
#include <helib/log.h>
#include <helib/keys.h>
#include <chrono>
#include <cstdlib>
#include <future>

using std::chrono::steady_clock;
using std::chrono::duration;
using std::chrono::duration_cast;

#ifdef HELIB_DEBUG

#include <helib/debugging.h>

namespace helib {

long printFlag = FLAG_PRINT_VEC;
/************************ Some local functions ***********************/
/*********************************************************************/

static void checkCriticalValue(const std::vector<NTL::ZZX>& zzParts,
                               const DoubleCRT& sKey,
                               const RecryptData& rcData,
                               long q);

static void checkRecryptBounds(const std::vector<NTL::ZZX>& zzParts,
                               const DoubleCRT& sKey,
                               const Context& context,
                               long q);

static void checkRecryptBounds_v(const std::vector<NTL::ZZX>& v,
                                 const DoubleCRT& sKey,
                                 const Context& context,
                                 long q);
} // namespace helib

#endif // HELIB_DEBUG

namespace helib {

// Return in poly a polynomial with X^i encoded in all the slots
static void x2iInSlots(NTL::ZZX& poly,
                       long i,
                       std::vector<NTL::ZZX>& xVec,
                       const EncryptedArray& ea)
{
  xVec.resize(ea.size());
  NTL::ZZX x2i = NTL::ZZX(i, 1);
  for (long j = 0; j < (long)xVec.size(); j++)
    xVec[j] = x2i;
  ea.encode(poly, xVec);
}

// Make every entry of vec divisible by p2e by adding/subtracting q, while
// keeping the added multiples small.  Specifically, for q = 1 mod p2e and any
// integer z can be made divisible by p2e via z' = z + v*q, with |v| <= p2e/2.

static void newMakeDivisible(NTL::ZZX& poly,
                             long p2e,
                             long q,
                             const Context& context,
                             NTL::ZZX& vpoly)
{
  if (p2e == 1) {
    vpoly = 0;
    return;
  }

  assertTrue<InvalidArgument>(q > 0l, "q must be positive");
  assertTrue<InvalidArgument>(p2e > 0l, "p2e must be positive");

  assertEq<InvalidArgument>(q % p2e, 1l, "q must equal 1 modulo p2e");

  long p = context.getP();

  const RecryptData& rcData = context.getRcData();
  const PowerfulDCRT& p2d_conv = *rcData.p2dConv;

  NTL::Vec<NTL::ZZ> pwrfl;
  p2d_conv.ZZXtoPowerful(pwrfl, poly);

#ifdef HELIB_DEBUG
  NTL::Vec<NTL::ZZ> vvec(NTL::INIT_SIZE, pwrfl.length());
#endif

  for (long i : range(pwrfl.length())) {
    NTL::ZZ& z = pwrfl[i];
    long v;

    // What to add to z to make it divisible by p2e?
    long zMod = rem(z, p2e); // zMod is in [0,p2e-1]
    // NOTE: this makes sure we get a truly balanced remainder
    if (zMod > p2e / 2 || (p == 2 && zMod == p2e / 2 && NTL::RandomBnd(2))) {
      // randomize so that v has expected value 0
      zMod = p2e - zMod;
    } else {
      // need to add a negative number
      zMod = -zMod;
    }
    v = zMod;
    z += NTL::to_ZZ(q) * v; // make z divisible by p2e

    if (rem(z, p2e) != 0) { // sanity check
      std::cerr << "**error: original z[" << i
                << "]=" << (z - (NTL::to_ZZ(q) * v)) << std::dec
                << ", p^e=" << p2e << std::endl;
      std::cerr << "z' = z + " << v << "*q = " << z << std::endl;
      exit(1);
    }

#ifdef HELIB_DEBUG
    vvec[i] = v;
#endif
  }

  p2d_conv.powerfulToZZX(poly, pwrfl);

#ifdef HELIB_DEBUG
  p2d_conv.powerfulToZZX(vpoly, vvec);
#endif
}

/**
 * store the scaled-up poly inside zzParts
 * store the overflow part inside I_part
 */
static void multAndGetOverflowPart(std::vector<NTL::ZZX>& zzParts,
                                   std::vector<NTL::ZZX>& I_part,
                                   const Context& context)
{
  assertTrue(zzParts.size() == 2, "expecting a two-part ctxt");
  long aux = context.getAux();
  // long eNew = context.getENew();
  // long p = context.getP();
  long qKS = context.ithPrime(context.getIndexQks());

  const RecryptData& rcData = context.getRcData();
  const PowerfulDCRT& p2d_conv = *rcData.p2dConv;

  NTL::Vec<NTL::ZZ> pwfl, pwfl_mod;
  I_part.resize(2);

  long halfQKS = (qKS + 1) >> 1;
  long halfAux = (aux + 1) >> 1;

  for (size_t i = 0; i < zzParts.size(); i++) {
    p2d_conv.ZZXtoPowerful(pwfl, zzParts[i]);
    long length = pwfl.length();
    pwfl_mod.SetLength(length);
    vecRed(pwfl, pwfl, qKS, false);
    for (int j = 0; j < length; j++) {
      pwfl[j] *= aux;
      // divide pwfl[i] by qKS, store the quotient in pwfl[i]
      // the remainder in pwfl_mod[i]
      pwfl_mod[j] = NTL::DivRem(pwfl[j], pwfl[j], qKS);
      if (pwfl_mod[j] >= halfQKS) { // reduce pwfl_mod mod qKS
        pwfl_mod[j] -= qKS;
        pwfl[j] += 1;
      }
      // aux*a = t*q + [aux*a]_q, t in pwfl and [*]_q in pwfl_mod
      // we actually want -t instead of t
      pwfl[j] = -pwfl[j];
      if (pwfl[j] >= halfAux) // reduce pwfl mod aux
        pwfl[j] -= aux;
    }
    // now pwfl_mod stores the multiplied & qKS-reduced powerful rep
    // pwfl stores the aux-reduced overflow part
    p2d_conv.powerfulToZZX(zzParts[i], pwfl_mod);
    p2d_conv.powerfulToZZX(I_part[i], pwfl);
    // NOTE: actually we should mod-up on the powerful basis
    //  i.e., mod-up pwfl_mod instead of zzParts
    //  but since the modulus is huge, this is still ok
  }
  // Ctxt will be constructed in PubKey::reCrypt, not here
}

/*********************************************************************/
/*********************************************************************/

/**
 * Summary of Appendix A from https://ia.cr/2014/873 (version from 2019):
 * Assume that we already chosen e, e' and t.
 *
 * Based in this analysis, we need
 *    (1) (f*p^{e'} + 2*p^r+2))*B <= p^e/2
 * where B is a certain high-probability bound and f is a certain
 * fudge factor.
 *
 **/

// the routine compute_fudge is used to correct for the fact that
// the v-coeffs are not quite uniform

static double compute_fudge(long p2ePrime, long p2e)
{
  double eps = 0;

  if (p2ePrime > 1) {

    if (p2ePrime % 2 == 0) {
      eps = 1 / fsquare(p2ePrime);

      // The exact variance in this case is at most the variance
      // of a random variable that is distributed over
      //    -N..+N
      // where N = 2^{e'}/2.
      // Each endpoint occurs with probability 1/(4*N),
      // and the remaining values each occur with the same probability
      // 1/(2*N)

      // This variance is exactly computed as
      //    (N^2)/3 + 1/6 = ((N^2)/3)*(1 + 1/(2*N^2)), where N = 2^{e'}/2
      // So the std dev is at most
      //    N/sqrt(3)*(1 + 1/(4*N^2))

    } else {
      eps = 1 / double(p2e);

      // We are computing X + Y mod p^{e'}, where
      // X and Y are independent.
      // Y is uniformly distributed over
      //    -floor(p^{r}/2)..floor(p^{r}/2)
      // X is distributed over
      //    -floor(p^e/2)-1..floor(p^e/2)+1,
      // where each endpoint occurs with probability 1 / (2*(p^e+1)),
      // and the remaining p^e values are equally likely

      // The variance in this case is bounded by
      //   (N^2)/3*(1-eps) + (N^2)*eps = (N^2)/3*(1+2*eps),
      //       where N = p^{e'}/2 and eps < 1/p^e
      // So the std dev is bounded by
      //    N/sqrt(3)*sqrt(1+2*eps) <= N/sqrt(3)*(1+eps)
    }
  }

  return 1 + eps;
}

double RecryptData::setAE(long& e, long& ePrime, const Context& context)
{
  double coeff_bound = context.boundForRecryption();
  // coeff_bound is ultimately a high prob bound on |w0+w1*s|,
  // the coeffs of w0, w1 are chosen uniformly on [-1/2,1/2]

  long p = context.getP();
  long p2r = context.getAlMod().getPPowR();
  long r = context.getAlMod().getR();
  long frstTerm = 2 * p2r + 2;

  long e_bnd = 0;
  long p2e_bnd = 1;
  while (p2e_bnd <= ((1L << 30) - 2) / p) { // NOTE: this avoids overflow
    e_bnd++;
    p2e_bnd *= p;
  }
  // e_bnd is largest e such that p^e+1 < 2^30

  // Start with the smallest e s.t. p^e/2 >= frstTerm*coeff_bound
  ePrime = 0;
  e = r + 1;
  while (e <= e_bnd && NTL::power_long(p, e) < frstTerm * coeff_bound * 2)
    e++;

  //  if (e > e_bnd) Error("setAE: cannot find suitable e");
  assertFalse<RuntimeError>(e > e_bnd, "setAE: cannot find suitable e");

  // long ePrimeTry = r+1;
  long ePrimeTry = 1;

  while (ePrimeTry <= e_bnd) {
    long p2ePrimeTry = NTL::power_long(p, ePrimeTry);
    // long eTry = ePrimeTry+1;
    long eTry = std::max(r + 1, ePrimeTry + 1);
    while (eTry <= e_bnd && eTry - ePrimeTry < e - ePrime) {
      long p2eTry = NTL::power_long(p, eTry);
      double fudge = compute_fudge(p2ePrimeTry, p2eTry);
      if (p2eTry >= (p2ePrimeTry * fudge + frstTerm) * coeff_bound * 2)
        break;

      eTry++;
    }

    if (eTry <= e_bnd && eTry - ePrimeTry < e - ePrime) {
      e = eTry;
      ePrime = ePrimeTry;
    }

    ePrimeTry++;
  }
  // set minCapacity
  double minCapacity = log(HELIB_MIN_CAP_FRAC * p2r * coeff_bound /
                           double(NTL::power_long(p, e) + 1) /
                           context.getZMStar().getNormBnd()) /
                       log(2.0);

#ifdef HELIB_DEBUG
  std::cerr << "RecryptData::setAE(): e=" << e << ", e'=" << ePrime
            << std::endl;
#endif
  return minCapacity;
}

double RecryptData::setEncapAE(long& e,
                               long& ePrime,
                               long& qks,
                               long& R,
                               long nDgtsBTS,
                               const Context& context)
{
  double coeff_bound = context.boundForRecryption();
  // coeff_bound is ultimately a high prob bound on |w0+w1*s|,
  // the coeffs of w0, w1 are chosen uniformly on [-1/2,1/2]

  long p = context.getP();
  long p2r = context.getAlMod().getPPowR();
  long r = context.getAlMod().getR();
  long frstTerm = 2 * p2r + 2;

  long e_bnd = 0;
  long p2e_bnd = 1;

  // while (p2e_bnd <= ((1L << 30) - 2) / p) { // NOTE: this avoids overflow
  //   e_bnd++;
  //   p2e_bnd *= p;
  // }
  while (p2e_bnd <=
         ((1L << NTL_SP_NBITS) - 2) / p) { // NOTE: this avoids overflow
    e_bnd++;
    p2e_bnd *= p;
  }
  // e_bnd is largest e such that p^e+1 < 2^30

  // Start with the smallest e s.t. p^e/2 >= frstTerm*coeff_bound
  ePrime = 0;
  e = r + 1;
  while (e <= e_bnd && NTL::power_long(p, e) < frstTerm * coeff_bound * 2)
    e++;

  //  if (e > e_bnd) Error("setAE: cannot find suitable e");
  assertFalse<RuntimeError>(e > e_bnd, "setAE: cannot find suitable e");

  // long ePrimeTry = r+1;
  long ePrimeTry = 1;

  while (ePrimeTry <= e_bnd) {
    long p2ePrimeTry = NTL::power_long(p, ePrimeTry);
    // long eTry = ePrimeTry+1;
    long eTry = std::max(r + 1, ePrimeTry + 1);
    while (eTry <= e_bnd && eTry - ePrimeTry < e - ePrime) {
      long p2eTry = NTL::power_long(p, eTry);
      double fudge = compute_fudge(p2ePrimeTry, p2eTry);
      if (p2eTry >= (p2ePrimeTry * fudge + frstTerm) * coeff_bound * 2)
        break;

      eTry++;
    }

    if (eTry <= e_bnd && eTry - ePrimeTry < e - ePrime) {
      e = eTry;
      ePrime = ePrimeTry;
    }

    ePrimeTry++;
  }

// #ifdef HELIB_DEBUG
  std::cerr << "RecryptData::setEncapAE(): e=" << e << ", e'=" << ePrime
            << std::endl;
// #endif

  // the logic for setting e and ePrime is the same,
  // the only difference is that we set qKS and R here
  // TODO: avoid copying the code of setAE
  double scale = context.getScale();
  long m = context.getZMStar().getM();
  long phim = context.getZMStar().getPhiM();
  long skHwt = context.getHwt();
  double mfac =
      context.getZMStar()
          .getNormBnd(); // Dm, the ratio between pwfl and canonical bound
  long q = NTL::power_long(p, e) + 1;
  double stdev = NTL::conv<double>(context.getStdev());

  // the "scaled noise" should be bounded by 2/3* p^r*(B*+0.5) (powerful)
  // still assuming the input ctxt before key-switching has a canonical noise
  // beta then after KS, noise is 2*beta*q/qKS (canonical) i.e.,
  // 2*beta*q/qKS*mfac < 2/3 * p^r(B*+0.5) or qKS > 2*beta*q*mfac / (2/3 *
  // p^r*(B*+0.5))
  double beta =
      p2r * (scale * sqrt(double(phim) / 12.0) *
                 (sqrt(double(skHwt) * log(double(phim))) + 1) +
             0.5); // beta is the canonical bound for the message after ms
  qks = ceil(2 * beta * q * mfac / (HELIB_MIN_CAP_FRAC * p2r * coeff_bound));
  printf("log2(qks) >= %f", log(double(qks)) / log(2.0));
  // now find the smallest qks with an m-th order primitive root of unity.
  // we only test for qks == 1 mod m
  long tmp = qks % m;
  if (tmp <= 1)
    qks += 1 - tmp;
  else
    qks += m + 1 - (qks % m);
  while (!NTL::ProbPrime(qks, 60))
    qks += m;
  printf(", final log2(qks) = %f\n", log(double(qks)) / log(2.0));
  // set the minimum capacity
  double minCapacity = log(qks / beta) / log(2.0);

  // finally, we choose R >= alpha / beta, where alpha is the bound on the
  // ks-added noise. for power-of-2 m, we have:
  // alpha = qKS^(1/c)*c*phim*sqrt(ln(phim))*P*sigma*k/sqrt(12)
  double alpha = pow(double(qks), 1.0 / nDgtsBTS) * phim *
                 sqrt(log(double(phim)) / 12.0) * p2r * scale * nDgtsBTS *
                 stdev;
  if (NTL::NumTwos(NTL::ZZ(m)) == 0) // alpha is larger for non-power-of-2 m
    alpha *= double(m) / sqrt(double(phim));
  R = ceil(alpha / beta);
  printf("log2(R) >= %f", log(double(R)) / log(2.0));
  // again, search for the next prime >= R with m-th root
  tmp = R % m;
  if (tmp <= 1)
    R += 1 - tmp;
  else
    R += m + 1 - tmp;
  while (!NTL::ProbPrime(R, 60))
    R += m;
  printf(", final log2(R) = %f\n", log(double(R)) / log(2.0));

  return minCapacity;
}

double RecryptData::setNewBootAE(long& eNew,
                                 long& t,
                                 long& aux,
                                 long& qks,
                                 long& R,
                                 long& numExtract,
                                 long nDgtsBTS,
                                 const Context& context)
{
  long p = context.getP();
  long p2r = context.getAlMod().getPPowR();
  long r = context.getAlMod().getR();
  double scale = context.getScale();
  long skHwt = context.getHwt();
  long encapSkHwt = context.getEncapHwt();
  long m = context.getZMStar().getM();
  long phim = context.getZMStar().getPhiM();
  double mfac =
      context.getZMStar()
          .getNormBnd(); // Dm, the ratio between pwfl and canonical bound
  long nfactors = context.getZMStar().getNFactors();
  double stdev = NTL::conv<double>(context.getStdev());

  // first,  we check the validity of t and set eNew & aux based on t
  double I_bound =
      context.boundForRecryption(); // B is the powerful bound for the overflow
                                    // part & ms-added noise (using the
                                    // encapsulated sparse key)
                                    // this implicitly use btsScale
  long I_range = 2 * long(ceil(I_bound)) + 1;
  // sanity check, no more than 2 digits to be extracted
  // TODO: maybe consider more than 2 digits later?
  assertTrue(I_range < p * p,
             "new bootstrap procedure supports extraction of at most 2 digits");
  if (t < 0) { // use non-power-of-p aux
    assertTrue(I_range * I_range < p,
               "non-power-of-p aux: p too small to hold the overflow part");
    aux = 2 * ceil(I_bound) + 1;
    if (NTL::GCD(aux, p) > 1)
      aux += 1;
    if (const char* explicitAux = std::getenv("HELIB_EXPLICIT_AUX")) {
      char* end = nullptr;
      long parsedAux = std::strtol(explicitAux, &end, 10);
      assertTrue(end != explicitAux && *end == '\0' && parsedAux > 1,
                 "invalid HELIB_EXPLICIT_AUX");
      assertTrue(NTL::GCD(parsedAux, p) == 1,
                 "HELIB_EXPLICIT_AUX must be coprime to p");
      aux = parsedAux;
      std::cout << "using HELIB_EXPLICIT_AUX=" << aux << "\n";
    }
    eNew = r;
    numExtract = 1;
  } else if (t == 0) { // deduce the value of t
    // we require p^t > I_range
    long tmp = 1;
    while (tmp <= I_range) {
      tmp *= p;
      t++;
    }
    eNew = r + t;
    aux = NTL::power_long(p, t);
    numExtract = t;
  } else { // check the validity of t
    assertTrue(
        NTL::power_long(p, t) > I_range,
        "preset power-of-p aux: p^t too small to hold the overflow part");
    eNew = r + t;
    aux = NTL::power_long(p, t);
    numExtract = ceil(log(I_range) / log(p) + 0.00001);
  }
  printf("\nt = %ld\n", t);
  printf("bound on I = %f\n", I_bound);
  printf("log2(aux) = %f\n", log(aux) / log(2));

  // second, we choose q
  // beta = P*(k*sqrt(phim/12)*(1 + sqrt(h*ln(phim))) + 0.5),
  // the canonical bound on the modulus switching noise
  double beta =
      p2r * (scale * sqrt(double(phim) / 12.0) *
                 (sqrt(double(skHwt) * log(double(phim))) + 1) +
             0.5); // beta is the canonical bound for the message after ms
  printf("log2(beta) = %f\n", log(beta) / log(2.0));
  // this is the pwfl bound on the modulus switching noise
  // 0.5 + k*2^(nfactors/2)/sqrt(12)*sqrt(encap_h)*sqrt(phim/m)
  double ms_bound =
      0.5 + scale *
                sqrt(double(phim) / double(m) * double(encapSkHwt) *
                     double(1L << nfactors) / 3.0) *
                0.5;
  // can2pwfl(2 * beta) + ms_pwfl_bound < q / (2 * aux)
  qks = ceil(aux * 2 * (mfac * 2 * beta + ms_bound * p2r + 2));
  printf("log2(qks) >= %f", log(double(qks)) / log(2.0));
  // now find the smallest qks with an m-th order primitive root of unity.
  // we only test for qks == 1 mod m
  long tmp = qks % m;
  if (tmp <= 1)
    qks += 1 - tmp;
  else
    qks += m + 1 - (qks % m);
  while (!NTL::ProbPrime(qks, 60))
    qks += m;
  printf(", final log2(qks) = %f\n", log(double(qks)) / log(2.0));
  // set min capacity
  double minCapacity = log(qks / beta) / log(2.0);

  // finally, we choose R >= alpha / beta, where alpha is the bound on the
  // ks-added noise. for power-of-2 m, we have:
  // alpha = qKS^(1/c)*c*phim*sqrt(ln(phim))*P*sigma*k/sqrt(12)
  double alpha = pow(double(qks), 1.0 / nDgtsBTS) * phim *
                 sqrt(log(double(phim)) / 12.0) * p2r * scale * nDgtsBTS *
                 stdev;
  if (NTL::NumTwos(NTL::ZZ(m)) == 0) // alpha is larger for non-power-of-2 m
    alpha *= double(m) / sqrt(double(phim));
  R = ceil(alpha / beta);
  printf("log2(R) >= %f", log(double(R)) / log(2.0));
  // again, search for the next prime >= R with m-th root
  tmp = R % m;
  if (tmp <= 1)
    R += 1 - tmp;
  else
    R += m + 1 - tmp;
  while (!NTL::ProbPrime(R, 60))
    R += m;
  printf(", final log2(R) = %f\n", log(double(R)) / log(2.0));

  return minCapacity;
}

bool RecryptData::operator==(const RecryptData& other) const
{
  if (mvec != other.mvec)
    return false;

  if (skHwt != other.skHwt)
    return false;

  return true;
}

// The main method
// TODO: collapsed FFT for power-of-2 full-SIMD parameters
void RecryptData::init(const Context& context,
                       const NTL::Vec<long>& mvec_,
                       bool enableThick,
                       bool build_cache_,
                       bool minimal)
{
  if (alMod != nullptr) { // were we called for a second time?
    std::cerr << "@Warning: multiple calls to RecryptData::init\n";
    return;
  }

  // sanity check
  assertEq(computeProd(mvec_),
           context.getM(),
           "Cyclotomic polynomial mismatch");

  // Record the arguments to this function
  mvec = mvec_;
  build_cache = build_cache_;
  alsoThick = enableThick;

  bool mvec_ok = true;
  for (long i : range(mvec.length())) {
    NTL::Vec<NTL::Pair<long, long>> factors;
    factorize(factors, mvec[i]);
    if (factors.length() > 1)
      mvec_ok = false;
  }

  if (!mvec_ok) {
    Warning("prime power factorization recommended for bootstrapping");
  }

  skHwt = context.getEncapHwt();
  e = context.getE();
  ePrime = context.getEPrime();
  // new bts
  eNew = context.getENew();
  t = context.getT();
  newBtsFlag = context.getNewBTSFlag();

  long r = context.getAlMod().getR();

  // First part of Bootstrapping works wrt plaintext space p^{r'}
  long new_hensel_lifting =
      newBtsFlag
          ? eNew
          : (e - ePrime +
             r); // NOTE: is it necessary to create a new AlMod when eNew == r?
  alMod =
      std::make_shared<PAlgebraMod>(context.getZMStar(), new_hensel_lifting);
  ea = std::make_shared<EncryptedArray>(context, *alMod);
  // Polynomial defaults to F0, PAlgebraMod explicitly given

  p2dConv = std::make_shared<PowerfulDCRT>(context, mvec);

  if (!enableThick)
    return;

  // Initialize the linear polynomial for unpacking the slots
  NTL::zz_pBak bak;
  bak.save();
  ea->getAlMod().restoreContext();
  long nslots = ea->size();
  long d = ea->getDegree();

  const NTL::Mat<NTL::zz_p>& CBi =
      ea->getDerived(PA_zz_p()).getNormalBasisMatrixInverse();

  std::vector<NTL::ZZX> LM;
  LM.resize(d);
  for (long i = 0; i < d; i++) // prepare the linear polynomial
    LM[i] = rep(CBi[i][0]);

  std::vector<NTL::ZZX> C;
  ea->buildLinPolyCoeffs(C, LM); // "build" the linear polynomial

  unpackSlotEncoding.resize(d); // encode the coefficients

  for (long j = 0; j < d; j++) {
    std::vector<NTL::ZZX> v(nslots);
    for (long k = 0; k < nslots; k++)
      v[k] = C[j];
    ea->encode(unpackSlotEncoding[j], v);
  }
  firstMap = std::make_shared<EvalMap>(*ea, minimal, mvec, true, build_cache);
  // for the new BTS, the secondMap is still performed modulo
  // p^new_hensel_lifting
  secondMap = std::make_shared<EvalMap>(newBtsFlag ? *ea : context.getEA(),
                                        minimal,
                                        mvec,
                                        false,
                                        build_cache);
}

/********************************************************************/
/********************************************************************/

// Extract digits from fully packed slots
// this is pre-declaration
BootBench extractDigitsPacked(Ctxt& ctxt,
                         long botHigh,
                         long r,
                         long ePrime,
                         const std::vector<NTL::ZZX>& unpackSlotEncoding);

// Extract digits from unpacked slots
void extractDigitsThin(Ctxt& ctxt, long botHigh, long r, long ePrime, bool thinRefine=false);

// #define HELIB_DEBUG // XXX

// bootstrap a ciphertext to reduce noise
BootBench PubKey::reCrypt(Ctxt& ctxt) const
{
  BootBench benchmarker;
  auto time_boot_start = steady_clock::now();
  HELIB_TIMER_START;

  // Some sanity checks for dummy ciphertext
  long ptxtSpace = ctxt.getPtxtSpace();
  if (ctxt.isEmpty())
    return benchmarker;
  if (ctxt.parts.size() == 1 && ctxt.parts[0].skHandle.isOne()) {
    // Dummy encryption, just ensure that it is reduced mod p
    NTL::ZZX poly = to_ZZX(ctxt.parts[0]);
    for (long i = 0; i < poly.rep.length(); i++)
      poly[i] = NTL::to_ZZ(rem(poly[i], ptxtSpace));
    poly.normalize();
    ctxt.DummyEncrypt(poly);
    return benchmarker;
  }

  // check that we have bootstrapping data
  assertTrue(recryptKeyID >= 0l, "No bootstrapping data");

  long p = getContext().getP();
  long r = getContext().getAlMod().getR();
  long p2r = getContext().getAlMod().getPPowR();

  long intFactor = ctxt.intFactor;

  // the bootstrapping key is encrypted relative to plaintext space p^{e-e'+r}.
  const RecryptData& rcData = getContext().getRcData();
  long e = rcData.e;
  long ePrime = rcData.ePrime;
  long p2ePrime = NTL::power_long(p, ePrime);
  long q = NTL::power_long(p, e) + 1;
  // new bts
  // TODO: remove the duplicate members in rcData to avoid possible mistakes?
  long eNew = rcData.eNew;
  long t = rcData.t;
  long aux = context.getAux();
  long p2eNew = NTL::power_long(p, eNew);
  bool newBtsFlag = context.getNewBTSFlag();
  bool newKSFlag = context.getNewKSFlag();
  if (!newBtsFlag)
    assertTrue(e >= r, "rcData.e must be at least alMod.r");
  if (newBtsFlag)
    q = context.ithPrime(context.getIndexQks());

#ifdef HELIB_DEBUG
  std::cerr << "reCrypt: p=" << p << ", r=" << r << ", e=" << e
            << " ePrime=" << ePrime << ", q=" << q << std::endl;
  CheckCtxt(ctxt, "init");
#endif

  // can only bootstrap ciphertext with plaintext-space dividing p^r
  assertEq(p2r % ptxtSpace, 0l, "ptxtSpace must divide p^r when bootstrapping");

#ifdef HELIB_DEBUG
  NTL::ZZX poly_input;
  dbgKey->Decrypt(poly_input, ctxt);
#endif

  ctxt.dropSmallAndSpecialPrimes();

#ifdef HELIB_DEBUG
  CheckCtxt(ctxt, "after drop small and special primes");
#endif

  HELIB_NTIMER_START(AAA_preProcess);

  // Make sure that this ciphertext is in canonical form
  if (!ctxt.inCanonicalForm())
    ctxt.reLinearize();

  // Mod-switch down if needed
  IndexSet s = ctxt.getPrimeSet() / context.getSpecialPrimes();
  assertTrue(s <= context.getCtxtPrimes(), "prime set is messed up");
  if (newKSFlag) { // since the decomposition is independent of the RNS
                   // representation, we do not need to keep 3 primes
    s.clear();
    s.insert(context.getIndexQks());
    ctxt.bringToSet(s);
    // TODO: check the error bound of ctxt and warn on large noise?
    //  the budget in ctxt should be larger than the budget in the
    //  assumed-post-ms-pre-ks-ctxt (whose noise bound is beta) by one or two
    //  bits
  } else {
    if (s.card() > 3) { // leave only first three ciphertext primes
      long first = s.first();
      IndexSet s3(first, first + 2);
      s.retain(s3);
    }
    ctxt.modDownToSet(s);
  }

#ifdef HELIB_DEBUG
  CheckCtxt(ctxt, "after mod down to KS modulus");

  NTL::ZZX poly_before;
  dbgKey->Decrypt(poly_before, ctxt);
  assertEq(poly_input, poly_before, "something wrong with mod-down");
#endif

  // key-switch to the bootstrapping key
  // > this is problematic... because the qKS will be dropped by
  // dropSmallAndSpecialPrimes
  // also, the decomposition algorithm needs to change
  ctxt.reLinearize(recryptKeyID, newKSFlag);

#ifdef HELIB_DEBUG
  CheckCtxt(ctxt, "after key switching");

  NTL::ZZX poly_after;
  dbgKey->Decrypt(poly_after, ctxt);
  assertEq(poly_before, poly_after, "something wrong with key switching");
#endif

  // "raw mod-switch" to the bootstrapping modulus q=p^e+1.
  std::vector<NTL::ZZX> zzParts; // the mod-switched parts, in ZZX format

  double mfac = ctxt.getContext().getZMStar().getNormBnd();
  // NOTE: rawModSwitch assumes the current ctxt modulus is coprime to q
  // however, this is not the case for the new bootstrap,
  // where ctxt modulus is qKS * R, while q = qKS
  double noise_est;
  if (newBtsFlag)
    noise_est = ctxt.rawModSwitchNew(zzParts, q) * mfac; // q = qKS
  else
    noise_est = ctxt.rawModSwitch(zzParts, q) * mfac; // q = p^e+1
  // noise_est is an upper bound on the L-infty norm of the scaled noise
  // in the pwrfl basis
  long phim = context.getZMStar().getPhiM();
  double beta =
      p2r * (context.getScale() * sqrt(double(phim) / 12.0) *
                 (sqrt(double(context.getHwt()) * log(double(phim))) + 1) +
             0.5);
  double noise_bnd = newBtsFlag ? 2 * beta * mfac
                                : HELIB_MIN_CAP_FRAC * p2r *
                                      ctxt.getContext().boundForRecryption();
  // noise_bnd is the bound assumed in selecting the parameters
  double noise_rat = noise_est / noise_bnd;

  HELIB_STATS_UPDATE("raw-mod-switch-noise", noise_rat);

  if (noise_rat > 1) {
    // TODO: Turn the following preprocessor logics into a warnOrThrow function
    std::string message =
        "rawModSwitch scaled noise exceeds bound: " + std::to_string(noise_rat);
#ifdef HELIB_DEBUG
    Warning(message);
#else
    throw LogicError(message);
#endif
  }

  assertEq(zzParts.size(),
           (std::size_t)2,
           "Exactly 2 parts required for mod-switching in thin bootstrapping");

#ifdef HELIB_DEBUG
  const PowerfulDCRT& p2d_conv = *rcData.p2dConv;
  if (dbgKey && !newBtsFlag) {
    checkRecryptBounds(zzParts, dbgKey->getRecryptKey(), ctxt.getContext(), q);
  }
  // after raw mod-swtich, poly_after has a stored intFactor of
  // `intFactor`, and an implied intFactor of `intFactor * [q]_p2r`
  NTL::ZZX poly_after_rawmod;
  NTL::vec_ZZ pwfl_after_rawmod;
  if (newBtsFlag) { // NOTE: this is pure debug
    rawDecrypt(poly_after_rawmod, zzParts, dbgKey->getRecryptKey());
    p2d_conv.ZZXtoPowerful(pwfl_after_rawmod, poly_after_rawmod);
    vecRed(pwfl_after_rawmod, pwfl_after_rawmod, q, false);
    p2d_conv.powerfulToZZX(poly_after_rawmod, pwfl_after_rawmod);
    // and may output incorrect results for out-of-range inputs
    long impliedIntFactor = NTL::MulMod(intFactor, q % p2r, p2r);
    long invImplied = NTL::InvMod(impliedIntFactor, p2r);
    // set abs=True to be consistent with the Decryption function
    poly_after_rawmod = helib::MulMod(poly_after_rawmod, invImplied, p2r, true);
    assertEq(poly_after,
             poly_after_rawmod,
             "something wrong with raw modswitch");
  }
#endif

  std::vector<NTL::ZZX> v;
  v.resize(2);

  // this ctxt will be
  //  (a) modulus switched down to the ctxtPrimes (b) key-switched to
  //  the normal sk (c) substracted by the digit extraction results
  Ctxt modUpCtxt(ctxt.getPubKey(), p2eNew);
  std::vector<NTL::ZZX> I_part;

  // the debug polys
  NTL::ZZX poly_before_mod;
  NTL::vec_ZZ pwfl_before_mod;
  NTL::ZZX poly_after_mod, I_poly;
  NTL::vec_ZZ pwfl_after_mod, pwfl_diff;
  NTL::ZZX modUpPoly;
  NTL::vec_ZZ modUpPwfl, I_part_pwfl;
  if (!newBtsFlag)
    // Add multiples of q to make the zzParts divisible by p^{e'}
    for (long i : range(2)) {
      // make divisible by p^{e'}
      newMakeDivisible(zzParts[i], p2ePrime, q, ctxt.getContext(), v[i]);
    }
  else {
    // XXX: modding up or down & relin should not affect the stored intFactor...
    assertTrue(intFactor == ctxt.intFactor, "why the intFactors do not match?");

#ifdef HELIB_DEBUG
    rawDecrypt(poly_before_mod, zzParts, dbgKey->getRecryptKey());
    p2d_conv.ZZXtoPowerful(pwfl_before_mod, poly_before_mod);
    vecRed(pwfl_before_mod, pwfl_before_mod, q, false);
    double largePwflBefore = NTL::conv<double>(largestCoeff(pwfl_before_mod));
    printf("pwfl capacity before = %f\n",
           log(double(q) / 2.0 / largePwflBefore / double(aux)) / log(2.0));
    assertTrue(largePwflBefore * aux < double(q) / 2.0,
               "not enough capacity before scaleUp");
#endif

    // > (1) mult by aux,
    // (2) extract the overflow part, and store them as ZZXs -> feed into digit
    // extraction (3) mod-up the ctxt after mult to the highest modulus (with an
    // extra prime modulus) note that the mod-up ctxt decrypts w.r.t. encapSk
    multAndGetOverflowPart(zzParts, I_part, ctxt.getContext());

#ifdef HELIB_DEBUG
    rawDecrypt(poly_after_mod, zzParts, dbgKey->getRecryptKey());
    // check the overflow part
    rawDecrypt(I_poly, I_part, dbgKey->getRecryptKey());
    p2d_conv.ZZXtoPowerful(I_part_pwfl, I_poly);
    vecRed(I_part_pwfl, I_part_pwfl, aux, false);
    p2d_conv.powerfulToZZX(I_poly, I_part_pwfl);

    p2d_conv.ZZXtoPowerful(pwfl_diff, poly_after_mod - I_poly * q);
    // vecRed(pwfl_diff, pwfl_diff, NTL::ZZ(aux) * NTL::ZZ(q), false);
    for (long i = 0; i < phim; i++) {
      assertTrue(bool(NTL::abs(pwfl_diff[i]) <= (q / 2)),
                 "overflow part wrong?");
    }
    // assertTrue(NTL::to_long(largestCoeff(pwfl_diff)) <= (q / 2),
    //            "overflow part wrong?");
    // check the relationship between pwfl_before_mod and pwfl_after_mod
    // now check if the bound on I is valid
    p2d_conv.ZZXtoPowerful(pwfl_after_mod, poly_after_mod);
    NTL::xdouble actualBoundI =
        NTL::conv<NTL::xdouble>(largestCoeff(pwfl_after_mod)) / NTL::xdouble(q);
    std::cerr << "actual bound on I is " << actualBoundI << "\n";
    assertTrue(
        bool(NTL::floor(actualBoundI) <= ceil(context.boundForRecryption())),
        "bound on I exceeded");
    vecRed(pwfl_after_mod, pwfl_after_mod, q, false);
    if (t < 0) {
      for (long i = 0; i < phim; i++) {
        long before = NTL::rem(pwfl_before_mod[i], p2r);
        before = NTL::MulMod(before, aux % p2r, p2r);
        long after = NTL::rem(pwfl_after_mod[i], p2r);
        assertEq(before,
                 after,
                 "something wrong with upscaling, non-power-of-p aux");
      }
    } else {
      for (long i = 0; i < phim; i++) {
        long before = NTL::rem(pwfl_before_mod[i], p2r);
        long after = NTL::rem(pwfl_after_mod[i], p2eNew);
        assertEq(before * aux,
                 after,
                 "something wrong with upscaling, power-of-p aux");
      }
    }
#endif

    // now handle the mod-up ctxt
    modUpCtxt.ptxtSpace = p2eNew;
    modUpCtxt.primeSet = context.getCtxtPrimes() | context.getModUpPrimes();
    long QmodP = 1;
    // be careful with single-precision modular arithmetic...
    for (auto i : modUpCtxt.primeSet)
      QmodP = NTL::MulMod(QmodP, context.ithPrime(i) % p2eNew, p2eNew);
    // we pretend the IMPLIED INTFACTOR of m* mod P is 1
    // (i.e., zzParts after raw-ms decrypts to m*)
    // which is actually intFactor * [q]_{p^r}
    // we want the mult-by-aux and modded-up ctxt to have
    // an implied intfactor of aux (decrypts to m*+...)
    if (t < 0) // gcd(aux, p) = 1, stored intFactor = aux * [Q^-1]_P
      modUpCtxt.intFactor =
          NTL::MulMod(aux, NTL::InvMod(QmodP, p2eNew), p2eNew);
    else // aux = p^t, stored intFactor = [Q^-1]_P
      modUpCtxt.intFactor = NTL::InvMod(QmodP, p2eNew);
    modUpCtxt.noiseBound = context.boundForRecryption() * q;
    for (int i = 0; i < 2; i++)
      modUpCtxt.addPart(DoubleCRT(zzParts[i], context, modUpCtxt.primeSet),
                        ctxt.parts[i].skHandle);
    // discard the mod-up primes
    modUpCtxt.modDownToSet(context.getCtxtPrimes());

    // switch to the normal sk (keyID = 0)
    modUpCtxt.reLinearize();

#ifdef HELIB_DEBUG
    dbgKey->Decrypt(modUpPoly, modUpCtxt);
    p2d_conv.ZZXtoPowerful(modUpPwfl, modUpPoly);
    vecRed(modUpPwfl, modUpPwfl, p2eNew, false);
    p2d_conv.ZZXtoPowerful(I_part_pwfl, I_poly);
    vecRed(I_part_pwfl, I_part_pwfl, aux, false);

    p2d_conv.ZZXtoPowerful(pwfl_after_mod, poly_after_mod);
    vecRed(pwfl_after_mod, pwfl_after_mod, q, false);

    // for aux = p^t
    // the message in zzParts befure "multAndGetOverflowPart" is m' mod p2r
    // after that, the scaled message becomes p^t*m' mod p2eNew
    // by modding up, the message becomes p^t*m' + [q]_p2eNew*I
    // setting the modUpCtxt.intFactor = [Q^-1]_p2eNew ensures
    // modUpCtxt decrypts to p^t*m' + [q]_p2eNew*I
    for (long i = 0; i < phim; i++) {
      if (t > 0) {
        assertEq(NTL::rem(modUpPwfl[i] - I_part_pwfl[i] * (q % p2eNew) -
                              pwfl_after_mod[i],
                          p2eNew),
                 0L,
                 "something wrong with assembly...? power-of-p aux");
      } else {
        // aux is not a power of p
        // old message: m' mod p2r
        // scaled message: aux*m' mod p2r (= pwfl_after_mod)
        // modding up: aux*m' + [q]_p2r*I
        // by setting stored intFactor to aux*[Q]_p2r:
        //  m' + [q*aux^-1]_p2r*I
        assertEq(NTL::rem(aux * modUpPwfl[i] - I_part_pwfl[i] * (q % p2r) -
                              pwfl_after_mod[i],
                          p2r),
                 0L,
                 "something wrong with assembly...? non-power-of-p aux");
      }
    }
#endif
  }

#ifdef HELIB_DEBUG
  if (dbgKey) {
    if (!newBtsFlag) {
      checkRecryptBounds_v(v, dbgKey->getRecryptKey(), ctxt.getContext(), q);
      checkCriticalValue(zzParts,
                         dbgKey->getRecryptKey(),
                         ctxt.getContext().getRcData(),
                         q);
    }
  }
#endif
  if (!newBtsFlag)
    for (long i : range(zzParts.size())) {
      zzParts[i] /= p2ePrime; // divide by p^{e'}
    }

  // NOTE: here we lose the intFactor associated with ctxt.
  // We will restore it below.
  ctxt = recryptEkey;

  // XXX: debug variables
  NTL::ZZX I_poly_new;
  NTL::vec_ZZ I_poly_new_pwfl;
  if (!newBtsFlag) {
    ctxt.multByConstant(zzParts[1]);
    ctxt.addConstant(zzParts[0]);
    benchmarker.bits_after_inner_prod = ctxt.capacity();
  } else {
    // this is ok, ctxt still decrypts to I + aux*I'
    ctxt.multByConstant(I_part[1]);
    ctxt.addConstant(I_part[0]);
    benchmarker.bits_after_inner_prod = ctxt.capacity();
#ifdef HELIB_DEBUG
    // NOTE-debug: check ctxt here, see if it matches I_part
    dbgKey->Decrypt(I_poly_new, ctxt);
    p2d_conv.ZZXtoPowerful(I_poly_new_pwfl, I_poly_new);
    vecRed(I_poly_new_pwfl, I_poly_new_pwfl, p2eNew, false);
    for (long i = 0; i < phim; i++) {
      assertTrue(bool(NTL::rem(I_poly_new_pwfl[i] - I_part_pwfl[i], aux) == 0),
                 "something wrong with linear dec");
    }
#endif
  }
#ifdef HELIB_DEBUG
  CheckCtxt(ctxt, "after preProcess");
#endif
  HELIB_NTIMER_STOP(AAA_preProcess);

#ifdef HELIB_DEBUG
  NTL::ZZX before_map;
  NTL::vec_ZZ before_map_pwfl;
  dbgKey->Decrypt(before_map, ctxt);
  p2d_conv.ZZXtoPowerful(before_map_pwfl, before_map);
  vecRed(before_map_pwfl, before_map_pwfl, p2eNew, false);
#endif

  // Move the powerful-basis coefficients to the plaintext slots
  HELIB_NTIMER_START(AAA_LinearTransform1);
  auto bits_down_linear1 = ctxt.capacity();
  auto time_linear1_start = steady_clock::now();
  ctxt.getContext().getRcData().firstMap->apply(ctxt);
  auto time_linear1_end = steady_clock::now();
  bits_down_linear1 -= ctxt.capacity();
  HELIB_NTIMER_STOP(AAA_LinearTransform1);

#ifdef HELIB_DEBUG
  CheckCtxt(ctxt, "after LinearTransform1");
#endif

  // Extract the digits e-e'+r-1,...,e-e' (from fully packed slots)
  HELIB_NTIMER_START(AAA_extractDigitsPacked);
  // > well, this needs to be heavily modified...
  // the case of new BTS is handled within the function, not here
  // auto bits_down_extract = ctxt.capacity();
  // auto time_extract_start = steady_clock::now();
  auto subBench = extractDigitsPacked(ctxt,
                      e - ePrime,
                      r,
                      ePrime,
                      context.getRcData().unpackSlotEncoding);
  // auto time_extract_end = steady_clock::now();
  // bits_down_extract -= ctxt.capacity();
  HELIB_NTIMER_STOP(AAA_extractDigitsPacked);

#ifdef HELIB_DEBUG
  CheckCtxt(ctxt, "after extractDigitsPacked");
#endif

  // Move the slots back to powerful-basis coefficients
  HELIB_NTIMER_START(AAA_LinearTransform2);
  auto bits_down_linear2 = ctxt.capacity();
  auto time_linear2_start = steady_clock::now();
  ctxt.getContext().getRcData().secondMap->apply(ctxt);
  auto time_linear2_end = steady_clock::now();
  bits_down_linear2 -= ctxt.capacity();
  HELIB_NTIMER_STOP(AAA_LinearTransform2);

#ifdef HELIB_DEBUG
  CheckCtxt(ctxt, "after linearTransform2");
  NTL::ZZX after_map;
  NTL::vec_ZZ after_map_pwfl;
  if (newBtsFlag) {
    dbgKey->Decrypt(after_map, ctxt);
    p2d_conv.ZZXtoPowerful(after_map_pwfl, after_map);
    vecRed(after_map_pwfl, after_map_pwfl, p2eNew, false);
    for (long i = 0; i < phim; i++) {
      assertEq(NTL::rem(before_map_pwfl[i] - after_map_pwfl[i], aux),
               0L,
               "ext step congruence failed");
      assertEq(after_map_pwfl[i], I_part_pwfl[i], "ext step reduction failed");
    }
  }
#endif

  // restore intFactor
  if (newBtsFlag) {
#ifdef HELIB_DEBUG
    { // XXX: debug, perform raw assemly
      // input message
      NTL::vec_ZZ input_pwfl;
      p2d_conv.ZZXtoPowerful(input_pwfl, poly_input);
      vecRed(input_pwfl, input_pwfl, p2r, false);
      // modUp message = modUpPwfl
      // afterMap message = after_map_pwfl
      if (t > 0) {
        for (long i = 0; i < phim; i++) {
          // NTL::ZZ p2eNewZZ(p2eNew);
          long tmp = NTL::rem(after_map_pwfl[i], p2eNew);
          tmp = NTL::MulMod(q % p2eNew, tmp, p2eNew);
          tmp = NTL::SubMod(NTL::rem(modUpPwfl[i], p2eNew), tmp, p2eNew);
          assertEq(tmp % aux, 0L, "not divisible by aux?");
          tmp /= aux;
          tmp %= p2r;
          long factorDiv = intFactor;
          factorDiv = NTL::MulMod(factorDiv % p2r, q % p2r, p2r);
          assertEq(tmp,
                   NTL::MulMod(NTL::rem(input_pwfl[i], p2r), factorDiv, p2r),
                   "not matching, p^t");
        }
      } else {
        assertEq(p2eNew, p2r, "??????");
        for (long i = 0; i < phim; i++) {
          long tmp = NTL::rem(after_map_pwfl[i], p2r);
          long mul_fac = NTL::MulMod(q % p2r, NTL::InvMod(aux % p2r, p2r), p2r);
          tmp = NTL::MulMod(tmp, mul_fac, p2r);
          tmp = NTL::SubMod(NTL::rem(modUpPwfl[i], p2r), tmp, p2r);

          long factorDiv = intFactor;
          factorDiv = NTL::MulMod(factorDiv % p2r, q % p2r, p2r);
          assertEq(tmp,
                   NTL::MulMod(NTL::rem(input_pwfl[i], p2r), factorDiv, p2r),
                   "not matching, delta_0");
        }
      }
    }
#endif
    // > assemble the ctxts and adjust the intFactors
    // after modding-up, the message encrypted becomes
    //  m + [q]_P*I, P = p^eNew, for power-of-p aux
    //  m + [q*aux^-1]_P*I, P = p^r, for non-power-of-p aux

    long adjust_factor = q % p2eNew;
    if (t < 0)
      adjust_factor = NTL::MulMod(q % p2r, NTL::InvMod(aux % p2r, p2r), p2r);
    ctxt *= adjust_factor;
    ctxt -= modUpCtxt;
    ctxt.negate(); // modUpCtxt - ctxt = Enc(m)
    // remove the extra p^t, if any
    for (long i = 0; i < t; i++)
      ctxt.divideByP();
    // now correct the intFactor
    // NOTE: for old bts, since q=p^e+1, [q]_p^r==1
    ctxt.intFactor =
        NTL::MulMod(ctxt.intFactor,
                    NTL::MulMod(intFactor, q % ptxtSpace, ptxtSpace),
                    ptxtSpace);
  } else if (intFactor != 1)
    ctxt.intFactor = NTL::MulMod(ctxt.intFactor, intFactor, ptxtSpace);
  auto time_boot_end = steady_clock::now();
  // benchmarker.bits_down_extract = bits_down_extract;
  benchmarker.bits_down_linear_1 = bits_down_linear1;
  benchmarker.bits_down_linear_2 = bits_down_linear2;
  benchmarker.bits_final = ctxt.capacity();
  benchmarker.time_linear_1 =
      duration_cast<duration<double>>(time_linear1_end - time_linear1_start)
          .count();
  benchmarker.time_linear_2 =
      duration_cast<duration<double>>(time_linear2_end - time_linear2_start)
          .count();
  // benchmarker.time_extract =
  //     duration_cast<duration<double>>(time_extract_end - time_extract_start)
  //         .count();
  benchmarker.time_total =
      duration_cast<duration<double>>(time_boot_end - time_boot_start).count();
  benchmarker += subBench;
  return benchmarker;
}

#ifdef HELIB_BOOT_THREADS

// Extract digits from fully packed slots, multithreaded version
// NOTE: this is invoked if HELIB_BOOT_THREADS is defined
BootBench extractDigitsPacked(Ctxt& ctxt,
                         long botHigh,
                         long r,
                         long ePrime,
                         const std::vector<NTL::ZZX>& unpackSlotEncoding)
{
  helib::BootBench benchmarker;
  HELIB_TIMER_START;

  // Step 1: unpack the slots of ctxt
  auto time_unpack_start = steady_clock::now();
  double cap_unpack_before = ctxt.capacity();
  HELIB_NTIMER_START(unpack);
  ctxt.cleanUp();

  // Apply the d automorphisms and store them in scratch area
  long d = ctxt.getContext().getOrdP();

  std::vector<Ctxt> unpacked(d, Ctxt(ZeroCtxtLike, ctxt));
  { // explicit scope to force all temporaries to be released
    std::vector<std::shared_ptr<DoubleCRT>> coeff_vector;
    std::vector<double> coeff_vector_sz;
    coeff_vector.resize(d);
    coeff_vector_sz.resize(d);

    HELIB_NTIMER_START(unpack1);
    for (long i = 0; i < d; i++) {
      coeff_vector[i] = std::make_shared<DoubleCRT>(unpackSlotEncoding[i],
                                                    ctxt.getContext(),
                                                    ctxt.getPrimeSet());
      coeff_vector_sz[i] = NTL::conv<double>(
          embeddingLargestCoeff(unpackSlotEncoding[i],
                                ctxt.getContext().getZMStar()));
    }
    HELIB_NTIMER_STOP(unpack1);

    HELIB_NTIMER_START(unpack2);
    std::vector<Ctxt> frob(d, Ctxt(ZeroCtxtLike, ctxt));

    NTL_EXEC_RANGE(d, first, last)
    // FIXME: implement using hoisting!
    for (long j = first; j < last; j++) { // process jth Frobenius
      frob[j] = ctxt;
      frob[j].frobeniusAutomorph(j);
      frob[j].cleanUp();
      // FIXME: not clear if we should call cleanUp here
    }
    NTL_EXEC_RANGE_END

    HELIB_NTIMER_STOP(unpack2);

    HELIB_NTIMER_START(unpack3);
    Ctxt tmp1(ZeroCtxtLike, ctxt);
    for (long i = 0; i < d; i++) {
      for (long j = 0; j < d; j++) {
        tmp1 = frob[j];
        tmp1.multByConstant(*coeff_vector[mcMod(i + j, d)],
                            coeff_vector_sz[mcMod(i + j, d)]);
        unpacked[i] += tmp1;
      }
    }
    HELIB_NTIMER_STOP(unpack3);
  }
  HELIB_NTIMER_STOP(unpack);
  auto time_unpack_end = steady_clock::now();
  double cap_unpack_after = cap_unpack_before;
  for (long i = 0; i < d; i++)
    cap_unpack_after = std::min(cap_unpack_after, unpacked[i].capacity());

  // #ifdef HELIB_DEBUG
  //   CheckCtxt(unpacked[0], "after unpack");
  // #endif

  double cap_ext_before = cap_unpack_after;
  auto time_ext_start = steady_clock::now();
  NTL_EXEC_RANGE(d, first, last)
  for (long i = first; i < last; i++) {
    extractDigitsThin(unpacked[i], botHigh, r, ePrime);
  }
  NTL_EXEC_RANGE_END
  auto time_ext_end = steady_clock::now();
  double cap_ext_after = cap_ext_before;
  for(long i = 0; i < d; i++)
    cap_ext_after = std::min(cap_ext_after, unpacked[i].capacity());

  // #ifdef HELIB_DEBUG
  //  CheckCtxt(unpacked[0], "before repack");
  // #endif

  bool newBtsFlag = ctxt.getContext().getNewBTSFlag();

  // Step 3: re-pack the slots
  double cap_repack_before = cap_ext_after;
  auto time_repack_start = steady_clock::now();
  HELIB_NTIMER_START(repack);
  const EncryptedArray& ea2 = newBtsFlag ? *ctxt.getContext().getRcData().ea
                                         : ctxt.getContext().getEA();
  NTL::ZZX xInSlots;
  std::vector<NTL::ZZX> xVec(ea2.size());
  ctxt = unpacked[0];
  for (long i = 1; i < d; i++) {
    x2iInSlots(xInSlots, i, xVec, ea2);
    unpacked[i].multByConstant(xInSlots);
    ctxt += unpacked[i];
  }
  HELIB_NTIMER_STOP(repack);
  auto time_repack_end = steady_clock::now();
  double cap_repack_after = ctxt.capacity();
  // #ifdef HELIB_DEBUG
  //  CheckCtxt(ctxt, "after repack");
  // #endif
  benchmarker.time_linear_1 = duration_cast<duration<double>>(time_unpack_end - time_unpack_start).count();
  benchmarker.time_linear_2 = duration_cast<duration<double>>(time_repack_end - time_repack_start).count();
  benchmarker.time_extract = duration_cast<duration<double>>(time_ext_end - time_ext_start).count();
  benchmarker.bits_down_linear_1 = cap_unpack_before - cap_unpack_after;
  benchmarker.bits_down_linear_2 = cap_repack_before - cap_repack_after;
  benchmarker.bits_down_extract = cap_ext_before - cap_ext_after;
  return benchmarker;
}

#else

// Extract digits from fully packed slots
// NOTE: this is invoked when HELIB_BOOT_THREADS is undefined
// botHigh = e - e' is the number of digits to discard
// r is the number of digits to keep
BootBench extractDigitsPacked(Ctxt& ctxt,
                         long botHigh,
                         long r,
                         long ePrime,
                         const std::vector<NTL::ZZX>& unpackSlotEncoding)
{
  helib::BootBench benchmarker;
  HELIB_TIMER_START;

  // Step 1: unpack the slots of ctxt
  auto time_unpack_start = steady_clock::now();
  double cap_unpack_before = ctxt.capacity();
  HELIB_NTIMER_START(unpack);
  ctxt.cleanUp();

  // Apply the d automorphisms and store them in scratch area
  long d = ctxt.getContext().getOrdP();

  std::vector<Ctxt> unpacked(d, Ctxt(ZeroCtxtLike, ctxt));
  { // explicit scope to force all temporaries to be released
    std::vector<std::shared_ptr<DoubleCRT>> coeff_vector;
    std::vector<double> coeff_vector_sz;
    coeff_vector.resize(d);
    coeff_vector_sz.resize(d);
    for (long i = 0; i < d; i++) {
      coeff_vector[i] = std::make_shared<DoubleCRT>(unpackSlotEncoding[i],
                                                    ctxt.getContext(),
                                                    ctxt.getPrimeSet());
      coeff_vector_sz[i] = NTL::conv<double>(
          embeddingLargestCoeff(unpackSlotEncoding[i],
                                ctxt.getContext().getZMStar()));
    }

    Ctxt tmp1(ZeroCtxtLike, ctxt);
    Ctxt tmp2(ZeroCtxtLike, ctxt);

    // FIXME: implement using hoisting!
    for (long j = 0; j < d; j++) { // process jth Frobenius
      tmp1 = ctxt;
      tmp1.frobeniusAutomorph(j);
      tmp1.cleanUp();
      // FIXME: not clear if we should call cleanUp here

      for (long i = 0; i < d; i++) {
        tmp2 = tmp1;
        tmp2.multByConstant(*coeff_vector[mcMod(i + j, d)],
                            coeff_vector_sz[mcMod(i + j, d)]);
        unpacked[i] += tmp2;
      }
    }
  }
  HELIB_NTIMER_STOP(unpack);
  auto time_unpack_end = steady_clock::now();
  double cap_unpack_after = cap_unpack_before;
  for (long i = 0; i < d; i++)
    cap_unpack_after = std::min(cap_unpack_after, unpacked[i].capacity());

  // #ifdef HELIB_DEBUG
  //   CheckCtxt(unpacked[0], "after unpack");
  // #endif

  double cap_ext_before = cap_unpack_after;
  auto time_ext_start = steady_clock::now();
  for (long i = 0; i < (long)unpacked.size(); i++) {
    extractDigitsThin(unpacked[i], botHigh, r, ePrime);
  }
  auto time_ext_end = steady_clock::now();
  double cap_ext_after = cap_ext_before;
  for(long i = 0; i < d; i++)
    cap_ext_after = std::min(cap_ext_after, unpacked[i].capacity());

  // #ifdef HELIB_DEBUG
  //   CheckCtxt(unpacked[0], "before repack");
  // #endif

  bool newBtsFlag = ctxt.getContext().getNewBTSFlag();

  // Step 3: re-pack the slots
  double cap_repack_before = cap_ext_after;
  auto time_repack_start = steady_clock::now();
  HELIB_NTIMER_START(repack);
  // NOTE: for the new bootstrap,
  //  the repacking still works on ptxt modulus p2eNew, not p2r
  const EncryptedArray& ea2 = newBtsFlag ? *ctxt.getContext().getRcData().ea
                                         : ctxt.getContext().getEA();
  NTL::ZZX xInSlots;
  std::vector<NTL::ZZX> xVec(ea2.size());
  ctxt = unpacked[0];
  for (long i = 1; i < d; i++) {
    x2iInSlots(xInSlots, i, xVec, ea2);
    unpacked[i].multByConstant(xInSlots);
    ctxt += unpacked[i];
  }
  HELIB_NTIMER_STOP(repack);
  auto time_repack_end = steady_clock::now();
  double cap_repack_after = ctxt.capacity();

  benchmarker.time_linear_1 = duration_cast<duration<double>>(time_unpack_end - time_unpack_start).count();
  benchmarker.time_linear_2 = duration_cast<duration<double>>(time_repack_end - time_repack_start).count();
  benchmarker.time_extract = duration_cast<duration<double>>(time_ext_end - time_ext_start).count();
  benchmarker.bits_down_linear_1 = cap_unpack_before - cap_unpack_after;
  benchmarker.bits_down_linear_2 = cap_repack_before - cap_repack_after;
  benchmarker.bits_down_extract = cap_ext_before - cap_ext_after;
  return benchmarker;
}

#endif

// Use packed bootstrapping, so we can bootstrap all in just one go.
void packedRecrypt(const CtPtrs& cPtrs,
                   const std::vector<zzX>& unpackConsts,
                   const EncryptedArray& ea)
{
  PubKey& pKey = (PubKey&)cPtrs[0]->getPubKey();

  // Allocate temporary ciphertexts for the recryption
  int nPacked = divc(cPtrs.size(), ea.getDegree()); // ceil(totalNum/d)
  std::vector<Ctxt> cts(nPacked, Ctxt(pKey));

  repack(CtPtrs_vectorCt(cts), cPtrs, ea); // pack ciphertexts
  //  cout << "@"<< lsize(cts)<<std::flush;
  for (Ctxt& c : cts) {   // then recrypt them
    c.reducePtxtSpace(2); // we only have recryption data for binary ctxt
    pKey.reCrypt(c);
  }
  unpack(cPtrs, CtPtrs_vectorCt(cts), ea, unpackConsts);
}

// recrypt all ctxt at level < belowLvl
void packedRecrypt(const CtPtrs& array,
                   const std::vector<zzX>& unpackConsts,
                   const EncryptedArray& ea,
                   long belowLvl)
{
  std::vector<Ctxt*> v;
  for (long i = 0; i < array.size(); i++)
    if (array.isSet(i) && !array[i]->isEmpty() &&
        array[i]->bitCapacity() < belowLvl * (array[i]->getContext().BPL()))
      v.push_back(array[i]);
  packedRecrypt(CtPtrs_vectorPt(v), unpackConsts, ea);
}
void packedRecrypt(const CtPtrMat& m,
                   const std::vector<zzX>& unpackConsts,
                   const EncryptedArray& ea,
                   long belowLvl)
{
  std::vector<Ctxt*> v;
  for (long i = 0; i < m.size(); i++)
    for (long j = 0; j < m[i].size(); j++)
      if (m[i].isSet(j) && !m[i][j]->isEmpty() &&
          m[i][j]->bitCapacity() < belowLvl * (m[i][j]->getContext().BPL()))
        v.push_back(m[i][j]);
  packedRecrypt(CtPtrs_vectorPt(v), unpackConsts, ea);
}

//===================== Thin Bootstrapping stuff ==================

void ThinRecryptData::init(const Context& context,
                           const NTL::Vec<long>& mvec_,
                           bool alsoThick,
                           bool build_cache_,
                           bool minimal)
{
  RecryptData::init(context, mvec_, alsoThick, build_cache_, minimal);
  coeffToSlot =
      std::make_shared<ThinEvalMap>(*ea, minimal, mvec, true, build_cache);
  // for thin bootstrapping
  // slotToCoeff takes place before I-part-extraction
  // thus it has the ordinary plaintext modulus
  slotToCoeff = std::make_shared<ThinEvalMap>(context.getEA(),
                                              minimal,
                                              mvec,
                                              false,
                                              build_cache);
}

// Extract digits from thinly packed slots

long fhe_force_chen_han = 0;

static bool auxCombinedPreSubtractEnabled()
{
  const char* flag = std::getenv("HELIB_AUX_COMBINED_PRE_SUBTRACT");
  return flag != nullptr && std::atol(flag) != 0;
}

static bool auxPreExtractSubtractEnabled()
{
  const char* flag = std::getenv("HELIB_AUX_PRE_EXTRACT_SUBTRACT");
  return flag != nullptr && std::atol(flag) != 0;
}

static bool auxBatchedCoeffToSlotEnabled()
{
  const char* flag = std::getenv("HELIB_AUX_BATCHED_COEFF2SLOT");
  return flag != nullptr && std::atol(flag) != 0;
}

static bool auxParallelCoeffToSlotEnabled()
{
  const char* flag = std::getenv("HELIB_AUX_PARALLEL_COEFF2SLOT");
  return flag != nullptr && std::atol(flag) != 0;
}

static bool auxLcAksDelayedCoeffToSlotEnabled()
{
  const char* flag = std::getenv("HELIB_AUX_LC_AKS_DELAYED");
  if (flag != nullptr && std::atol(flag) != 0)
    return true;
  flag = std::getenv("HELIB_AUX_LC_AKS_BSGS");
  if (flag != nullptr && std::atol(flag) != 0)
    return true;
  flag = std::getenv("HELIB_AUX_LC_AKS_FUSED_DIGITS");
  const char* allow = std::getenv("HELIB_AUX_LC_AKS_ALLOW_UNVERIFIED");
  return flag != nullptr && std::atol(flag) != 0 && allow != nullptr &&
         std::atol(allow) != 0;
}

// botHigh = e - e' is the number of digits to discard
// r is the number of digits to keep
void extractDigitsThin(Ctxt& ctxt, long botHigh, long r, long ePrime, bool thinRefine)
{
  HELIB_TIMER_START;

  Ctxt unpacked(ctxt);
  unpacked.cleanUp();

  std::vector<Ctxt> scratch;

  long p = ctxt.getContext().getP();
  long p2r = NTL::power_long(p, r);
  long topHigh =
      botHigh + r - 1; // topHigh is the index of the highest digit to keep

  bool newBtsFlag = ctxt.getContext().getNewBTSFlag();

  // degree Chen/Han technique is p^{bot-1}(p-1)r
  // degree of basic technique is p^{bot-1}p^r,
  //     or p^{bot-1}p^{r-1} if p==2, r > 1, and bot+r > 2

  if (newBtsFlag) {
    newExtractDigits(scratch, unpacked);
    long scratch_len = scratch.size();
    if(thinRefine) {
      // just remove the lower digits
      for (long i = 0; i < scratch_len; i++) {
        unpacked -= scratch[i];
        unpacked.divideByP();
      }
      ctxt = unpacked;
    } else {
      // > get the lower digits modulo p^eNew
      // sum up the extracted digits from high to low
      unpacked = scratch[scratch_len - 1];
      for (long i = scratch_len - 2; i >= 0; i--) {
        unpacked.multByP();
        unpacked += scratch[i];
      }
      ctxt = unpacked;
    }
#ifdef HELIB_DEBUG
    // XXX: debug
    // long p2eNew = NTL::power_long(p, ctxt.getContext().getENew());
    // std::vector<NTL::ZZX> in_slots, after_slots0;
    // ctxt.getContext().getRcData().ea->decrypt(ctxt, *dbgKey, in_slots);

    // ctxt.getContext().getRcData().ea->decrypt(scratch[0],
    //                                           *dbgKey,
    //                                           after_slots0);

    // long nslots = in_slots.size();
    // for (long i = 0; i < nslots; i++) {
    //   long in_deg = NTL::deg(in_slots[i]);
    //   long out_deg = NTL::deg(after_slots0[i]);
    //   assertTrue((in_deg <= 0) && (out_deg <= 0),
    //              "slot values should be integers");
    //   if ((in_deg == -1 && out_deg != -1) || (in_deg != -1 && out_deg == -1))
    //     assertTrue(false, "nope");
    //   if (in_deg == -1 && out_deg == -1)
    //     continue;
    //   NTL::ZZ tmp = in_slots[i][0];
    //   // long tmp_bal = NTL::rem(tmp, p2eNew);
    //   // if(tmp_bal > p2eNew / 2)
    //   //   tmp_bal -= p2eNew;
    //   long expected = NTL::rem(tmp, p);
    //   expected = helib::balRem(expected, p);
    //   long got = NTL::to_long(after_slots0[i][0]);
    //   got = helib::balRem(got, p2eNew);
    //   assertEq(expected, got, "???");
    // }
    // XXX: end debug
#endif
    return;
  }

  // XXX: pre-computed in Context
  // bool use_chen_han = false;
  // if (r > 1) {
  //   double chen_han_cost = log(p - 1) + log(r);
  //   double basic_cost;
  //   if (p == 2 && r > 2 && botHigh + r > 2)
  //     basic_cost = (r - 1) * log(p);
  //   else
  //     basic_cost = r * log(p);

  //   // std::cerr << "*** basic: " << basic_cost << "\n";
  //   // std::cerr << "*** chen/han: " << chen_han_cost << "\n";

  //   double thresh = 1.5;
  //   if (p == 2)
  //     thresh = 1.75;
  //   // increasing thresh makes chen_han less likely to be chosen.
  //   // For p == 2, the basic algorithm is just squaring,
  //   // and so is a bit cheaper, so we raise thresh a bit.
  //   // This is all a bit heuristic.

  //   if (basic_cost > thresh * chen_han_cost)
  //     use_chen_han = true;
  // }

  // if (fhe_force_chen_han > 0)
  //   use_chen_han = true;
  // else if (fhe_force_chen_han < 0)
  //   use_chen_han = false;

  if (ctxt.getContext().getCH18Flag()) {
    // use Chen and Han technique

    extendExtractDigits(scratch, unpacked, botHigh, r);

#if 0
    for (long i: range(scratch.size())) {
      CheckCtxt(scratch[i], "**");
    }
#endif

    for (long j = 0; j < botHigh; j++) {
      unpacked -= scratch[j];
      unpacked.divideByP();
    }

    if (p == 2 && botHigh > 0) // For p==2, subtract also the previous bit
      unpacked += scratch[botHigh - 1];
    unpacked.negate();

    if (r > ePrime) { // Add in digits from the bottom part, if any
      long topLow = r - 1 - ePrime;
      Ctxt tmp = scratch[topLow];
      for (long j = topLow - 1; j >= 0; --j) {
        tmp.multByP();
        tmp += scratch[j];
      }
      if (ePrime > 0)
        tmp.multByP(ePrime); // multiply by p^e'
      unpacked += tmp;
    }
    unpacked.reducePtxtSpace(p2r); // Our plaintext space is now mod p^r

    ctxt = unpacked;
  } else {

    if (p == 2 && r > 2 && topHigh + 1 > 2)
      topHigh--; // For p==2 we sometime get a bit for free

    extractDigits(scratch, unpacked, topHigh + 1);

    // set unpacked = -\sum_{j=botHigh}^{topHigh} scratch[j] * p^{j-botHigh}
    if (topHigh >= LONG(scratch.size())) {
      topHigh = scratch.size() - 1;
      std::cerr << " @ suspect: not enough digits in extractDigitsPacked\n";
    }

    unpacked = scratch[topHigh];
    for (long j = topHigh - 1; j >= botHigh; --j) {
      unpacked.multByP();
      unpacked += scratch[j];
    }
    if (p == 2 && botHigh > 0) // For p==2, subtract also the previous bit
      unpacked += scratch[botHigh - 1];
    unpacked.negate();

    if (r > ePrime) { // Add in digits from the bottom part, if any
      long topLow = r - 1 - ePrime;
      Ctxt tmp = scratch[topLow];
      for (long j = topLow - 1; j >= 0; --j) {
        tmp.multByP();
        tmp += scratch[j];
      }
      if (ePrime > 0)
        tmp.multByP(ePrime); // multiply by p^e'
      unpacked += tmp;
    }
    unpacked.reducePtxtSpace(p2r); // Our plaintext space is now mod p^r
    ctxt = unpacked;
  }
}

// Hack to get at private fields of public key
struct PubKeyHack
{                         // The public key
  const Context& context; // The context

  //! @var Ctxt pubEncrKey
  //! The public encryption key is an encryption of 0,
  //! relative to the first secret key
  Ctxt pubEncrKey;

  std::vector<long> skHwts;            // The Hamming weight of the secret keys
  std::vector<KeySwitch> keySwitching; // The key-switching matrices

  // The keySwitchMap structure contains pointers to key-switching matrices
  // for re-linearizing automorphisms. The entry keySwitchMap[i][n] contains
  // the index j such that keySwitching[j] is the first matrix one needs to
  // use when re-linearizing s_i(X^n).
  std::vector<std::vector<long>> keySwitchMap;

  NTL::Vec<int> KS_strategy; // NTL Vec's support I/O, which is
                             // more convenient

  // bootstrapping data

  long recryptKeyID; // index of the bootstrapping key
  Ctxt recryptEkey;  // the key itself, encrypted under key #0
};

#ifdef HELIB_DEBUG
static NTL::vec_ZZ balrem_vec_ZZX(const std::vector<NTL::ZZX>& vec, long q)
{
  NTL::vec_ZZ ret_vec;
  ret_vec.SetLength(vec.size());
  for (size_t i = 0; i < vec.size(); i++) {
    if (NTL::deg(vec[i]) == -1)
      ret_vec[i] = 0;
    else
      ret_vec[i] = balRem(NTL::rem(vec[i][0], q), q);
  }
  return ret_vec;
}
#endif


// bootstrap a ciphertext to reduce noise
BootBench PubKey::thinReCrypt(Ctxt& ctxt) const
{
#ifdef HELIB_DEBUG
  // XXX: debug, check the bound (5) in HS21
  if(context.getNewKSFlag())
  {
    long phim = context.getPhiM();
    long q = 1 << 20;
    NTL::zz_p::init(q);
    NTL::ZZX dbg_b, dbg_a, dbg_m;
    NTL::vec_ZZ pwfl_b, pwfl_a, pwfl_m;
    NTL::zz_p tmp;
    pwfl_b.SetLength(phim);
    pwfl_a.SetLength(phim);
    long n_trials = 20;
    auto p2d_conv = context.getRcData().p2dConv;
    NTL::ZZ sum_var(0);
    for (long i = 0; i < n_trials; i++) {
      // random in [-q/2,q/2]
      for (long j = 0; j < phim; j++) {
        NTL::random(tmp);
        pwfl_a[j] = NTL::conv<NTL::ZZ>(tmp);
        NTL::random(tmp);
        pwfl_b[j] = NTL::conv<NTL::ZZ>(tmp);
      }
      vecRed(pwfl_a, pwfl_a, q, false);
      vecRed(pwfl_b, pwfl_b, q, false);
      // check bound
      p2d_conv->powerfulToZZX(dbg_b, pwfl_b);
      p2d_conv->powerfulToZZX(dbg_a, pwfl_a);
      rawDecrypt(dbg_m, {dbg_b, dbg_a}, dbgKey->getRecryptKey());
      p2d_conv->ZZXtoPowerful(pwfl_m, dbg_m);
      NTL::ZZ pwfl_bound = largestCoeff(pwfl_m);
      std::cout << "pwfl bound on I " << 
                NTL::conv<NTL::xdouble>(pwfl_bound) / NTL::conv<NTL::xdouble>(q)
                << "\n";
      NTL::ZZ local_var(0), local_mean(0);
      for (long j = 0; j < phim; j++) {
        local_var += pwfl_m[j] * pwfl_m[j];
        local_mean += pwfl_m[j];
      }
      // += phim^2*Var(Iq)
      sum_var += phim * local_var - local_mean * local_mean;
    }
    // now sum_var = n_trials * phim^2 * Var(Iq)
    // std(I) = sqrt(Var(I)) = sqrt(sum_var / (n_trials * phim^2 * q^2))
    std::cout << "standard deviation "
              << NTL::sqrt(NTL::conv<NTL::xdouble>(sum_var) / n_trials) /
                     (q * phim)
              << "\n";
  }
#endif

  BootBench benchmarker;
  auto time_boot_start = steady_clock::now();
  HELIB_TIMER_START;

  // Some sanity checks for dummy ciphertext
  long ptxtSpace = ctxt.getPtxtSpace();
  if (ctxt.isEmpty())
    return benchmarker;

  if (ctxt.parts.size() == 1 && ctxt.parts[0].skHandle.isOne()) {
    // Dummy encryption, just ensure that it is reduced mod p
    NTL::ZZX poly = to_ZZX(ctxt.parts[0]);
    for (long i = 0; i < poly.rep.length(); i++)
      poly[i] = NTL::to_ZZ(rem(poly[i], ptxtSpace));
    poly.normalize();
    ctxt.DummyEncrypt(poly);
    return benchmarker;
  }

  // check that we have bootstrapping data
  assertTrue(recryptKeyID >= 0l, "Bootstrapping data not present");

  long p = ctxt.getContext().getP();
  long r = ctxt.getContext().getAlMod().getR();
  long p2r = ctxt.getContext().getAlMod().getPPowR();

  long intFactor = ctxt.intFactor;

  const ThinRecryptData& trcData = ctxt.getContext().getRcData();
  auto p2d_conv = *trcData.p2dConv;

  // the bootstrapping key is encrypted relative to plaintext space p^{e-e'+r}.
  long e = trcData.e;
  long ePrime = trcData.ePrime;
  long p2ePrime = NTL::power_long(p, ePrime);
  long q = NTL::power_long(p, e) + 1;
  // new bts
  long phim = context.getZMStar().getPhiM();
#ifdef HELIB_DEBUG
  long d = context.getOrdP();
  long n_slots = phim / d;
#endif
  long eNew = trcData.eNew;
  long t = trcData.t;
  long aux = context.getAux();
  long p2eNew = NTL::power_long(p, eNew);
  bool newBtsFlag = context.getNewBTSFlag();
  bool newKSFlag = context.getNewKSFlag();
  if (!newBtsFlag)
    assertTrue(e >= r, "trcData.e must be at least alMod.r");
  if (newBtsFlag)
    q = context.ithPrime(context.getIndexQks());

  // can only bootstrap ciphertext with plaintext-space dividing p^r
  assertEq(p2r % ptxtSpace,
           0l,
           "ptxtSpace must divide p^r when thin bootstrapping");

#ifdef HELIB_DEBUG
  auto ea_ctxt = context.shareEA();
  auto ea_boot = context.getRcData().ea;
  std::vector<NTL::ZZX> slots_in;
  ea_ctxt->decrypt(ctxt, *dbgKey, slots_in);
  CheckCtxt(ctxt, "init");
#endif

  ctxt.dropSmallAndSpecialPrimes();

#define DROP_BEFORE_THIN_RECRYPT
#define THIN_RECRYPT_NLEVELS (3)
#ifdef DROP_BEFORE_THIN_RECRYPT
  // experimental code...we should drop down to a reasonably low level
  // before doing the first linear map.
  long first = context.getCtxtPrimes().first();
  long last = std::min(context.getCtxtPrimes().last(),
                       first + THIN_RECRYPT_NLEVELS - 1);
  ctxt.bringToSet(IndexSet(first, last));
#endif

#ifdef HELIB_DEBUG
  CheckCtxt(ctxt, "after mod down");
#endif

  // Move the slots to powerful-basis coefficients
  HELIB_NTIMER_START(AAA_slotToCoeff);
  auto bits_down_linear1 = ctxt.capacity();
  auto time_linear1_start = steady_clock::now();
  trcData.slotToCoeff->apply(ctxt);
  auto time_linear1_end = steady_clock::now();
  bits_down_linear1 -= ctxt.capacity();
  HELIB_NTIMER_STOP(AAA_slotToCoeff);
  // std::cout << "time for linear1 is " << 
  //   duration_cast<duration<double>>(time_linear1_end - time_linear1_start).count() << "\n";

#ifdef HELIB_DEBUG
  CheckCtxt(ctxt, "after slotToCoeff");
#endif

  HELIB_NTIMER_START(AAA_bootKeySwitch);

  // Make sure that this ciphertext is in canonical form
  if (!ctxt.inCanonicalForm())
    ctxt.reLinearize();
  
#ifdef HELIB_DEBUG
  // XXX: debug
  NTL::ZZX poly_s2c;
  NTL::vec_ZZ pwfl_s2c;
  dbgKey->Decrypt(poly_s2c, ctxt);
  p2d_conv.ZZXtoPowerful(pwfl_s2c, poly_s2c);
  // NOTE: the pwfl basis coeffs stores
  // slots[0][0], ..., slots[0][d-1]; slots[1][0], ...
  vecRed(pwfl_s2c, pwfl_s2c, p2r, true);
  // check the validity of s2c
  for(long i = 0; i < n_slots; i++) {
    NTL::ZZ val_in;
    if(NTL::deg(slots_in[i]) == -1)
      val_in = 0;
    else
      val_in = NTL::rem(slots_in[i][0], p2r);
    assertEq(val_in, pwfl_s2c[i * d], "something wrong with s2c?");
  }
  // try
  // for(long i = 0; i < n_slots; i++) {
  //   auto it = std::find(pwfl_s2c.begin(), pwfl_s2c.end(), NTL::ZZ(i+1));
  //   if(it != pwfl_s2c.end())
  //     printf("%ld found at %ld\n", i+1, it - pwfl_s2c.begin());
  // }
  // XXX: end debug
#endif

  // Mod-switch down if needed
  IndexSet s = ctxt.getPrimeSet() / context.getSpecialPrimes();
  assertTrue(s <= context.getCtxtPrimes(), "prime set is messed up");
  if (newKSFlag) {
    s.clear();
    s.insert(context.getIndexQks());
    ctxt.bringToSet(s);
  } else {
    if (s.card() > 3) { // leave only first three ciphertext primes
      long first = s.first();
      IndexSet s3(first, first + 2);
      s.retain(s3);
    }
    ctxt.modDownToSet(s);
  }


#ifdef HELIB_DEBUG
  CheckCtxt(ctxt, "after mod down to qKS");
  // XXX: debug
  NTL::ZZX poly_beforeKS;
  NTL::vec_ZZ pwfl_beforeKS;
  dbgKey->Decrypt(poly_beforeKS, ctxt);
  p2d_conv.ZZXtoPowerful(pwfl_beforeKS, poly_beforeKS);
  vecRed(pwfl_beforeKS, pwfl_beforeKS, p2r, true);
  for(long i = 0; i < n_slots; i++) {
    assertEq(pwfl_s2c[i * d], pwfl_beforeKS[i * d], "something wrong with mod down to qKS");
  }
  // XXX: end debug
#endif


  // key-switch to the bootstrapping key
  ctxt.reLinearize(recryptKeyID, newKSFlag);

#ifdef HELIB_DEBUG
  CheckCtxt(ctxt, "after key switching");
  // XXX: debug
  NTL::ZZX poly_beforeRaw;
  NTL::vec_ZZ pwfl_beforeRaw;
  dbgKey->Decrypt(poly_beforeRaw, ctxt);
  p2d_conv.ZZXtoPowerful(pwfl_beforeRaw, poly_beforeRaw);
  vecRed(pwfl_beforeRaw, pwfl_beforeRaw, p2r, true);
  for(long i = 0; i < n_slots; i++) {
    assertEq(pwfl_beforeKS[i * d], pwfl_beforeRaw[i * d], "something wrong with ks to btk");
  }
  // XXX: end debug
#endif


  // "raw mod-switch" to the bootstrapping mosulus q=p^e+1.
  std::vector<NTL::ZZX> zzParts; // the mod-switched parts, in ZZX format

  double mfac = ctxt.getContext().getZMStar().getNormBnd();
  double noise_est;
  if (newBtsFlag)
    noise_est = ctxt.rawModSwitchNew(zzParts, q) * mfac;
  else
    noise_est = ctxt.rawModSwitch(zzParts, q) * mfac;
#ifdef HELIB_DEBUG
  // XXX: debug, check the bound on I after native bts
  NTL::ZZX native_poly;
  NTL::vec_ZZ native_pwfl;
  rawDecrypt(native_poly, zzParts, dbgKey->getRecryptKey());
  p2d_conv.ZZXtoPowerful(native_pwfl, native_poly);
  auto native_pwfl_bound = NTL::conv<NTL::xdouble>(largestCoeff(native_pwfl)) /
                           NTL::conv<NTL::xdouble>(q);
  std::cout << "native pwfl bound is " << native_pwfl_bound << "\n";
  // XXX: end debug native
#endif

  // noise_est is an upper bound on the L-infty norm of the scaled noise
  // in the pwrfl basis
  double beta =
      p2r * (context.getScale() * sqrt(double(phim) / 12.0) *
                 (sqrt(double(context.getHwt()) * log(double(phim))) + 1) +
             0.5);
  double noise_bnd = newBtsFlag ? 2 * beta * mfac
                                : HELIB_MIN_CAP_FRAC * p2r *
                                      ctxt.getContext().boundForRecryption();
  // noise_bnd is the bound assumed in selecting the parameters
  double noise_rat = noise_est / noise_bnd;

  HELIB_STATS_UPDATE("raw-mod-switch-noise", noise_rat);

  if (noise_rat > 1) {
    // TODO: Turn the following preprocessor logics into a warnOrThrow function
    std::string message =
        "rawModSwitch scaled noise exceeds bound: " + std::to_string(noise_rat);
#ifdef HELIB_DEBUG
    Warning(message);
#else
    throw LogicError(message);
#endif
  }

  assertEq(zzParts.size(),
           (std::size_t)2,
           "Exactly 2 parts required for mod-switching in thin bootstrapping");

#ifdef HELIB_DEBUG
  if (dbgKey) {
    if (!newBtsFlag)
      checkRecryptBounds(zzParts,
                         dbgKey->getRecryptKey(),
                         ctxt.getContext(),
                         q);
  }
#endif

  std::vector<NTL::ZZX> v;
  v.resize(2);

  Ctxt modUpCtxt(ctxt.getPubKey(), p2eNew);
  std::vector<NTL::ZZX> I_part;
  bool auxCombinedPreSubtract =
      newBtsFlag && t < 0 && auxCombinedPreSubtractEnabled();
  bool auxPreExtractSubtract =
      newBtsFlag && t < 0 && !auxCombinedPreSubtract &&
      auxPreExtractSubtractEnabled();
  bool auxBatchedCoeffToSlot =
      newBtsFlag && t < 0 && !auxCombinedPreSubtract &&
      !auxPreExtractSubtract && auxBatchedCoeffToSlotEnabled();
  bool auxParallelCoeffToSlot =
      newBtsFlag && t < 0 && !auxCombinedPreSubtract &&
      !auxPreExtractSubtract && !auxBatchedCoeffToSlot &&
      auxParallelCoeffToSlotEnabled();
  bool auxLcAksDelayedCoeffToSlot =
      newBtsFlag && t < 0 && !auxCombinedPreSubtract &&
      !auxPreExtractSubtract && !auxBatchedCoeffToSlot &&
      !auxParallelCoeffToSlot && auxLcAksDelayedCoeffToSlotEnabled();

  // debug variables

  NTL::ZZX poly_before_mod;
  NTL::vec_ZZ pwfl_before_mod;
  NTL::ZZX poly_after_mod, I_poly, modUpPoly;
  NTL::vec_ZZ I_part_pwfl, pwfl_diff, pwfl_after_mod, modUpPwfl;
  std::vector<NTL::ZZX> slots_scaleUp;

  steady_clock::time_point time_linear2_start_second = steady_clock::now();
  steady_clock::time_point time_linear2_end_second = time_linear2_start_second;
  if (!newBtsFlag) {
    // Add multiples of q to make the zzParts divisible by p^{e'}
    for (long i : range(2)) {
      // make divisible by p^{e'}

      newMakeDivisible(zzParts[i], p2ePrime, q, ctxt.getContext(), v[i]);
    }
  } else {
#ifdef HELIB_DEBUG
    // check the capacity of zzParts
    rawDecrypt(poly_before_mod, zzParts, dbgKey->getRecryptKey());
    p2d_conv.ZZXtoPowerful(pwfl_before_mod, poly_before_mod);
    vecRed(pwfl_before_mod, pwfl_before_mod, q, false);
    NTL::ZZ maxPwfl = largestCoeff(pwfl_before_mod);
    double capacity =
        log(double(q) / NTL::to_long(maxPwfl) / 2.0 / aux) / log(2.0);
    printf("pwfl capacity of zzParts is %f\n", capacity);
    assertTrue(capacity > 0, "not enough capacity before scale");
    // check the correctness of raw MS
    long divFac = NTL::MulMod(q % p2r, intFactor, p2r);
    divFac = NTL::InvMod(divFac, p2r);
    for(long i = 0; i < n_slots; i++) {
      long before = NTL::to_long(pwfl_beforeRaw[i * d]);
      long after = NTL::rem(pwfl_before_mod[i * d], p2r);
      after = NTL::MulMod(after, divFac, p2r);
      assertEq(after, before, "something wrong with raw ms");
    }
#endif
    // NOTE: in thin bootstrapping, the Zpr values in slots are first
    //  move to the coeff domain by slotToCoeff
    //  then the multAndGetOverflowPart adds noise to the lower digits
    //  and the unoccupied coeffs
    //  a coeffToSlot map moves the noisy coeffs into the slots,
    //  then a trace-like map removes the values out of Z_p2eNew
    //  In other words, we have apply the second map to
    //  both modUpCtxt and Enc(I+aux*I', p2eNew)
    multAndGetOverflowPart(zzParts, I_part, ctxt.getContext());

#ifdef HELIB_DEBUG
    // XXX: debug
    rawDecrypt(poly_after_mod, zzParts, dbgKey->getRecryptKey());
    // check the overflow part
    rawDecrypt(I_poly, I_part, dbgKey->getRecryptKey());
    p2d_conv.ZZXtoPowerful(I_part_pwfl, I_poly);
    vecRed(I_part_pwfl, I_part_pwfl, aux, false);
    p2d_conv.powerfulToZZX(I_poly, I_part_pwfl);

    p2d_conv.ZZXtoPowerful(pwfl_diff, poly_after_mod - I_poly * q);
    // vecRed(pwfl_diff, pwfl_diff, NTL::ZZ(aux) * NTL::ZZ(q), false);
    for (long i = 0; i < phim; i++) {
      assertTrue(bool(NTL::abs(pwfl_diff[i]) <= (q / 2)),
                 "overflow part wrong?");
    }
    // assertTrue(NTL::to_long(largestCoeff(pwfl_diff)) <= (q / 2),
    //            "overflow part wrong?");
    // check the relationship between pwfl_before_mod and pwfl_after_mod
    // now check if the bound on I is valid
    p2d_conv.ZZXtoPowerful(pwfl_after_mod, poly_after_mod);
    NTL::xdouble actualBoundI =
        NTL::conv<NTL::xdouble>(largestCoeff(pwfl_after_mod)) / NTL::xdouble(q);
    std::cerr << "actual bound on I is " << actualBoundI << "\n";
    assertTrue(
        bool(NTL::floor(actualBoundI) <= ceil(context.boundForRecryption())),
        "bound on I exceeded");
    vecRed(pwfl_after_mod, pwfl_after_mod, q, false);
    if (t < 0) {
      for (long i = 0; i < phim; i++) {
        long before = NTL::rem(pwfl_before_mod[i], p2r);
        before = NTL::MulMod(before, aux % p2r, p2r);
        long after = NTL::rem(pwfl_after_mod[i], p2r);
        assertEq(before,
                 after,
                 "something wrong with upscaling, non-power-of-p aux");
      }
    } else {
      for (long i = 0; i < phim; i++) {
        long before = NTL::rem(pwfl_before_mod[i], p2r);
        long after = NTL::rem(pwfl_after_mod[i], p2eNew);
        assertEq(before * aux,
                 after,
                 "something wrong with upscaling, power-of-p aux");
      }
    }
// XXX: end debug
#endif

    // now handle the mod-up ctxt
    modUpCtxt.ptxtSpace = p2eNew;
    modUpCtxt.primeSet = context.getCtxtPrimes() | context.getModUpPrimes();
    long QmodP = 1;
    // be careful with single-precision modular arithmetic...
    for (auto i : modUpCtxt.primeSet)
      QmodP = NTL::MulMod(QmodP, context.ithPrime(i) % p2eNew, p2eNew);
    // we pretend the IMPLIED INTFACTOR of m* mod P is 1
    // (i.e., zzParts after raw-ms decrypts to m*)
    // which is actually intFactor * [q]_{p^r}
    // we want the mult-by-aux and modded-up ctxt to have
    // an implied intfactor of aux (decrypts to m*+...)
    if (t < 0) // gcd(aux, p) = 1, stored intFactor = aux * [Q^-1]_P
      modUpCtxt.intFactor =
          NTL::MulMod(aux, NTL::InvMod(QmodP, p2eNew), p2eNew);
    else // aux = p^t, stored intFactor = [Q^-1]_P
      modUpCtxt.intFactor = NTL::InvMod(QmodP, p2eNew);
    modUpCtxt.noiseBound = context.boundForRecryption() * q;
    for (int i = 0; i < 2; i++)
      modUpCtxt.addPart(DoubleCRT(zzParts[i], context, modUpCtxt.primeSet),
                        ctxt.parts[i].skHandle);
    // discard the mod-up primes
    modUpCtxt.modDownToSet(context.getCtxtPrimes());

    // Switch to the normal sk (keyID = 0).  In the experimental delayed
    // LC-AKS route we keep modUpCtxt under the sparse boot key through the
    // first coeffToSlot transform, then switch back immediately after that
    // transform.  This is a correctness probe for the sparse-to-dense
    // through-linear-transform direction inspired by LCR+AKS.
    if (!auxLcAksDelayedCoeffToSlot)
      modUpCtxt.reLinearize();

#ifdef HELIB_DEBUG
    dbgKey->Decrypt(modUpPoly, modUpCtxt);
    p2d_conv.ZZXtoPowerful(modUpPwfl, modUpPoly);
    vecRed(modUpPwfl, modUpPwfl, p2eNew, false);
    p2d_conv.ZZXtoPowerful(I_part_pwfl, I_poly);
    vecRed(I_part_pwfl, I_part_pwfl, aux, false);

    p2d_conv.ZZXtoPowerful(pwfl_after_mod, poly_after_mod);
    vecRed(pwfl_after_mod, pwfl_after_mod, q, false);

    // for aux = p^t
    // the message in zzParts befure "multAndGetOverflowPart" is m' mod p2r
    // after that, the scaled message becomes p^t*m' mod p2eNew
    // by modding up, the message becomes p^t*m' + [q]_p2eNew*I
    // setting the modUpCtxt.intFactor = [Q^-1]_p2eNew ensures
    // modUpCtxt decrypts to p^t*m' + [q]_p2eNew*I
    for (long i = 0; i < phim; i++) {
      if (t > 0) {
        assertEq(NTL::rem(modUpPwfl[i] - I_part_pwfl[i] * (q % p2eNew) -
                              pwfl_after_mod[i],
                          p2eNew),
                 0L,
                 "something wrong with assembly...? power-of-p aux");
      } else {
        // aux is not a power of p
        // old message: m' mod p2r
        // scaled message: aux*m' mod p2r (= pwfl_after_mod)
        // modding up: aux*m' + [q]_p2r*I
        // by setting stored intFactor to aux*[Q]_p2r:
        //  m' + [q*aux^-1]_p2r*I
        assertEq(NTL::rem(aux * modUpPwfl[i] - I_part_pwfl[i] * (q % p2r) -
                              pwfl_after_mod[i],
                          p2r),
                 0L,
                 "something wrong with assembly...? non-power-of-p aux");
      }
    }
#endif
    if (auxCombinedPreSubtract) {
      Ctxt uncleanI(recryptEkey);
      uncleanI.multByConstant(I_part[1]);
      uncleanI.addConstant(I_part[0]);
      long adjust_factor =
          NTL::MulMod(q % p2r, NTL::InvMod(aux % p2r, p2r), p2r);
      uncleanI *= adjust_factor;
      modUpCtxt -= uncleanI;
      ctxt = modUpCtxt;
      benchmarker.bits_after_inner_prod = ctxt.capacity();
    } else if (!auxPreExtractSubtract && !auxBatchedCoeffToSlot &&
               !auxParallelCoeffToSlot) {
      // NOTE: this map is exclusive to thin-bts
      time_linear2_start_second = steady_clock::now();
      trcData.coeffToSlot->apply(modUpCtxt);
      if (auxLcAksDelayedCoeffToSlot)
        modUpCtxt.reLinearize();
      time_linear2_end_second = steady_clock::now();
      // std::cout << "time for 2nd linear2 is " << 
      //   duration_cast<duration<double>>(time_linear2_end_second - time_linear2_start_second).count()
      //   << "\n";
    }
#ifdef HELIB_DEBUG
    ea_boot->decrypt(modUpCtxt, *dbgKey, slots_scaleUp);
#endif
  }

#ifdef HELIB_DEBUG
  if (dbgKey) {
    if (!newBtsFlag) {
      checkRecryptBounds_v(v, dbgKey->getRecryptKey(), ctxt.getContext(), q);
      checkCriticalValue(zzParts,
                         dbgKey->getRecryptKey(),
                         ctxt.getContext().getRcData(),
                         q);
    }
  }
#endif

  if (!newBtsFlag)
    for (long i : range(zzParts.size())) {
      zzParts[i] /= p2ePrime; // divide by p^{e'}
    }

  // NOTE: here we lose the intFactor associated with ctxt.
  // We will restore it below.
  if (!auxCombinedPreSubtract)
    ctxt = recryptEkey;

  NTL::ZZX I_poly_new;
  NTL::vec_ZZ I_poly_new_pwfl;
  if (!newBtsFlag) {
    ctxt.multByConstant(zzParts[1]);
    ctxt.addConstant(zzParts[0]);
    benchmarker.bits_after_inner_prod = ctxt.capacity();
  } else if (!auxCombinedPreSubtract) {
    ctxt.multByConstant(I_part[1]);
    ctxt.addConstant(I_part[0]);
    benchmarker.bits_after_inner_prod = ctxt.capacity();
#ifdef HELIB_DEBUG
    // NOTE-debug: check ctxt here, see if it matches I_part
    dbgKey->Decrypt(I_poly_new, ctxt);
    p2d_conv.ZZXtoPowerful(I_poly_new_pwfl, I_poly_new);
    vecRed(I_poly_new_pwfl, I_poly_new_pwfl, p2eNew, false);
    for (long i = 0; i < phim; i++) {
      assertTrue(bool(NTL::rem(I_poly_new_pwfl[i] - I_part_pwfl[i], aux) == 0),
                 "something wrong with linear dec");
    }
#endif
  }
#ifdef HELIB_DEBUG
  CheckCtxt(ctxt, "after bootKeySwitch");
#endif

  HELIB_NTIMER_STOP(AAA_bootKeySwitch);

  double bits_down_extract = 0.0;
  steady_clock::time_point time_extract_start = steady_clock::now();
  steady_clock::time_point time_extract_end = time_extract_start;

  if (auxPreExtractSubtract) {
    HELIB_NTIMER_START(AAA_extractDigitsThin);
    bits_down_extract = ctxt.capacity();
    time_extract_start = steady_clock::now();
    extractDigitsThin(ctxt, e - ePrime, r, ePrime);
    time_extract_end = steady_clock::now();
    bits_down_extract -= ctxt.capacity();
    HELIB_NTIMER_STOP(AAA_extractDigitsThin);

    long adjust_factor =
        NTL::MulMod(q % p2r, NTL::InvMod(aux % p2r, p2r), p2r);
    ctxt *= adjust_factor;
    ctxt -= modUpCtxt;
    ctxt.negate();
  }

  // Move the powerful-basis coefficients to the plaintext slots
  HELIB_NTIMER_START(AAA_coeffToSlot);
  auto bits_down_linear2 = ctxt.capacity();
  auto time_linear2_start = steady_clock::now();
  if (auxBatchedCoeffToSlot)
    trcData.coeffToSlot->apply(ctxt, modUpCtxt);
  else if (auxParallelCoeffToSlot) {
    auto modUpFuture = std::async(std::launch::async, [&]() {
      trcData.coeffToSlot->apply(modUpCtxt);
    });
    trcData.coeffToSlot->apply(ctxt);
    modUpFuture.get();
  }
  else
    trcData.coeffToSlot->apply(ctxt);
  auto time_linear2_end = steady_clock::now();
  bits_down_linear2 -= ctxt.capacity();
  HELIB_NTIMER_STOP(AAA_coeffToSlot);
  // std::cout << "time for linear2 is " << 

#ifdef HELIB_DEBUG
  std::vector<NTL::ZZX> slots_I;
  ea_boot->decrypt(ctxt, *dbgKey, slots_I);
  // XXX: debug, check if the slots match
  auto vals_in = balrem_vec_ZZX(slots_in, p2r);
  auto vals_scaleUp = balrem_vec_ZZX(slots_scaleUp, p2eNew);
  auto val_I = balrem_vec_ZZX(slots_I, p2eNew);
  if (t > 0)
    for (long i = 0; i < n_slots; i++) {
      long tmp = balRem(NTL::rem(val_I[i], aux), aux); // I without I*
      tmp = (tmp + p2eNew) % p2eNew;
      tmp = NTL::MulMod(q % p2eNew, tmp, p2eNew);
      tmp = NTL::SubMod(NTL::rem(vals_scaleUp[i], p2eNew), tmp, p2eNew);
      assertEq(tmp % aux, 0L, "not divisible by aux?");
      tmp /= aux;
      tmp %= p2r;
      long factorDiv = intFactor;
      factorDiv = NTL::MulMod(factorDiv % p2r, q % p2r, p2r);
      assertEq(tmp,
               NTL::MulMod(NTL::rem(vals_in[i], p2r), factorDiv, p2r),
               "not matching, p^t");
    }
    // TODO: t < 0
    // XXX: end debug
#endif

#ifdef HELIB_DEBUG
  CheckCtxt(ctxt, "after coeffToSlot");
#endif

  // Extract the digits e-e'+r-1,...,e-e' (from fully packed slots)
  if (!auxPreExtractSubtract) {
    HELIB_NTIMER_START(AAA_extractDigitsThin);
    bits_down_extract = ctxt.capacity();
    time_extract_start = steady_clock::now();
    if (!auxCombinedPreSubtract)
      extractDigitsThin(ctxt, e - ePrime, r, ePrime);
    time_extract_end = steady_clock::now();
    bits_down_extract -= ctxt.capacity();
    HELIB_NTIMER_STOP(AAA_extractDigitsThin);
  }

#ifdef HELIB_DEBUG
  if (newBtsFlag) {
    std::vector<NTL::ZZX> slots_I_out;
    ea_boot->decrypt(ctxt, *dbgKey, slots_I_out);
    auto vals_I_out = balrem_vec_ZZX(slots_I_out, p2eNew);
    for (long i = 0; i < n_slots; i++) {
      assertEq(balRem(NTL::rem(val_I[i], aux), aux),
               NTL::to_long(vals_I_out[i]),
               "something wrong with digit extraction");
    }
  }
#endif

#ifdef HELIB_DEBUG
  CheckCtxt(ctxt, "after extractDigitsThin");
#endif

  // restore intFactor
  if (newBtsFlag) {
    // > assemble the ctxts and adjust the intFactors
    // after modding-up, the message encrypted becomes
    //  m + [q]_P*I, P = p^eNew, for power-of-p aux
    //  m + [q*aux^-1]_P*I, P = p^r, for non-power-of-p aux
    if (!auxCombinedPreSubtract && !auxPreExtractSubtract) {
      long adjust_factor = q % p2eNew;
      if (t < 0)
        adjust_factor = NTL::MulMod(q % p2r, NTL::InvMod(aux % p2r, p2r), p2r);
      ctxt *= adjust_factor;
      ctxt -= modUpCtxt;
      ctxt.negate();               // modUpCtxt - ctxt = Enc(m)
      for (long i = 0; i < t; i++) // remove the extra p^t, if any
        ctxt.divideByP();
    }
    // now correct the intFactor
    // NOTE: for old bts, since q=p^e+1, [q]_p^r==1
    ctxt.intFactor =
        NTL::MulMod(ctxt.intFactor,
                    NTL::MulMod(intFactor, q % ptxtSpace, ptxtSpace),
                    ptxtSpace);
  } else if (intFactor != 1)
    ctxt.intFactor = NTL::MulMod(ctxt.intFactor, intFactor, ptxtSpace);
  auto time_boot_end = steady_clock::now();
  benchmarker.bits_down_extract = bits_down_extract;
  benchmarker.bits_down_linear_1 = bits_down_linear1;
  benchmarker.bits_down_linear_2 = bits_down_linear2;
  benchmarker.bits_final = ctxt.capacity();
  benchmarker.time_linear_1 =
      duration_cast<duration<double>>(time_linear1_end - time_linear1_start)
          .count();
  benchmarker.time_linear_2 =
      duration_cast<duration<double>>(time_linear2_end - time_linear2_start)
          .count();
  if (newBtsFlag)
    benchmarker.time_linear_2 += 
      duration_cast<duration<double>>(time_linear2_end_second - time_linear2_start_second).count();
  benchmarker.time_extract =
      duration_cast<duration<double>>(time_extract_end - time_extract_start)
          .count();
  benchmarker.time_total =
      duration_cast<duration<double>>(time_boot_end - time_boot_start).count();
  return benchmarker;
}


// bootstrap a ciphertext to reduce noise
BootBench PubKey::thinReCryptRefine(Ctxt& ctxt) const
{
#ifdef HELIB_DEBUG
  // XXX: debug, check the bound (5) in HS21
  if(context.getNewKSFlag())
  {
    long phim = context.getPhiM();
    long q = 1 << 20;
    NTL::zz_p::init(q);
    NTL::ZZX dbg_b, dbg_a, dbg_m;
    NTL::vec_ZZ pwfl_b, pwfl_a, pwfl_m;
    NTL::zz_p tmp;
    pwfl_b.SetLength(phim);
    pwfl_a.SetLength(phim);
    long n_trials = 20;
    auto p2d_conv = context.getRcData().p2dConv;
    NTL::ZZ sum_var(0);
    for (long i = 0; i < n_trials; i++) {
      // random in [-q/2,q/2]
      for (long j = 0; j < phim; j++) {
        NTL::random(tmp);
        pwfl_a[j] = NTL::conv<NTL::ZZ>(tmp);
        NTL::random(tmp);
        pwfl_b[j] = NTL::conv<NTL::ZZ>(tmp);
      }
      vecRed(pwfl_a, pwfl_a, q, false);
      vecRed(pwfl_b, pwfl_b, q, false);
      // check bound
      p2d_conv->powerfulToZZX(dbg_b, pwfl_b);
      p2d_conv->powerfulToZZX(dbg_a, pwfl_a);
      rawDecrypt(dbg_m, {dbg_b, dbg_a}, dbgKey->getRecryptKey());
      p2d_conv->ZZXtoPowerful(pwfl_m, dbg_m);
      NTL::ZZ pwfl_bound = largestCoeff(pwfl_m);
      std::cout << "pwfl bound on I " << 
                NTL::conv<NTL::xdouble>(pwfl_bound) / NTL::conv<NTL::xdouble>(q)
                << "\n";
      NTL::ZZ local_var(0), local_mean(0);
      for (long j = 0; j < phim; j++) {
        local_var += pwfl_m[j] * pwfl_m[j];
        local_mean += pwfl_m[j];
      }
      // += phim^2*Var(Iq)
      sum_var += phim * local_var - local_mean * local_mean;
    }
    // now sum_var = n_trials * phim^2 * Var(Iq)
    // std(I) = sqrt(Var(I)) = sqrt(sum_var / (n_trials * phim^2 * q^2))
    std::cout << "standard deviation "
              << NTL::sqrt(NTL::conv<NTL::xdouble>(sum_var) / n_trials) /
                     (q * phim)
              << "\n";
  }
#endif

  BootBench benchmarker;
  auto time_boot_start = steady_clock::now();
  HELIB_TIMER_START;

  // Some sanity checks for dummy ciphertext
  long ptxtSpace = ctxt.getPtxtSpace();
  if (ctxt.isEmpty())
    return benchmarker;

  if (ctxt.parts.size() == 1 && ctxt.parts[0].skHandle.isOne()) {
    // Dummy encryption, just ensure that it is reduced mod p
    NTL::ZZX poly = to_ZZX(ctxt.parts[0]);
    for (long i = 0; i < poly.rep.length(); i++)
      poly[i] = NTL::to_ZZ(rem(poly[i], ptxtSpace));
    poly.normalize();
    ctxt.DummyEncrypt(poly);
    return benchmarker;
  }

  // check that we have bootstrapping data
  assertTrue(recryptKeyID >= 0l, "Bootstrapping data not present");

  long p = ctxt.getContext().getP();
  long r = ctxt.getContext().getAlMod().getR();
  long p2r = ctxt.getContext().getAlMod().getPPowR();

  long intFactor = ctxt.intFactor;

  const ThinRecryptData& trcData = ctxt.getContext().getRcData();
  auto p2d_conv = *trcData.p2dConv;

  // the bootstrapping key is encrypted relative to plaintext space p^{e-e'+r}.
  long e = trcData.e;
  long ePrime = trcData.ePrime;
  long p2ePrime = NTL::power_long(p, ePrime);
  long q = NTL::power_long(p, e) + 1;
  // new bts
  long phim = context.getZMStar().getPhiM();
#ifdef HELIB_DEBUG
  long d = context.getOrdP();
  long n_slots = phim / d;
  long aux = context.getAux();
#endif
  long eNew = trcData.eNew;
  long t = trcData.t;
  assertTrue(t > 0, "thin bts without I applies only to Delta=p^t");
  long p2eNew = NTL::power_long(p, eNew);
  bool newBtsFlag = context.getNewBTSFlag();
  assertTrue(newBtsFlag, "thinRecryptRefine expects new bts flag");
  bool newKSFlag = context.getNewKSFlag();
  if (!newBtsFlag)
    assertTrue(e >= r, "trcData.e must be at least alMod.r");
  if (newBtsFlag)
    q = context.ithPrime(context.getIndexQks());

  // can only bootstrap ciphertext with plaintext-space dividing p^r
  assertEq(p2r % ptxtSpace,
           0l,
           "ptxtSpace must divide p^r when thin bootstrapping");

#ifdef HELIB_DEBUG
  auto ea_ctxt = context.shareEA();
  auto ea_boot = context.getRcData().ea;
  std::vector<NTL::ZZX> slots_in;
  ea_ctxt->decrypt(ctxt, *dbgKey, slots_in);
  CheckCtxt(ctxt, "init");
#endif

  ctxt.dropSmallAndSpecialPrimes();

#define DROP_BEFORE_THIN_RECRYPT
#define THIN_RECRYPT_NLEVELS (3)
#ifdef DROP_BEFORE_THIN_RECRYPT
  // experimental code...we should drop down to a reasonably low level
  // before doing the first linear map.
  long first = context.getCtxtPrimes().first();
  long last = std::min(context.getCtxtPrimes().last(),
                       first + THIN_RECRYPT_NLEVELS - 1);
  ctxt.bringToSet(IndexSet(first, last));
#endif

#ifdef HELIB_DEBUG
  CheckCtxt(ctxt, "after mod down");
#endif

  // Move the slots to powerful-basis coefficients
  HELIB_NTIMER_START(AAA_slotToCoeff);
  auto bits_down_linear1 = ctxt.capacity();
  auto time_linear1_start = steady_clock::now();
  trcData.slotToCoeff->apply(ctxt);
  auto time_linear1_end = steady_clock::now();
  bits_down_linear1 -= ctxt.capacity();
  HELIB_NTIMER_STOP(AAA_slotToCoeff);
  // std::cout << "time for linear1 is " << 
  //   duration_cast<duration<double>>(time_linear1_end - time_linear1_start).count() << "\n";

#ifdef HELIB_DEBUG
  CheckCtxt(ctxt, "after slotToCoeff");
#endif

  HELIB_NTIMER_START(AAA_bootKeySwitch);

  // Make sure that this ciphertext is in canonical form
  if (!ctxt.inCanonicalForm())
    ctxt.reLinearize();
  
#ifdef HELIB_DEBUG
  // XXX: debug
  NTL::ZZX poly_s2c;
  NTL::vec_ZZ pwfl_s2c;
  dbgKey->Decrypt(poly_s2c, ctxt);
  p2d_conv.ZZXtoPowerful(pwfl_s2c, poly_s2c);
  // NOTE: the pwfl basis coeffs stores
  // slots[0][0], ..., slots[0][d-1]; slots[1][0], ...
  vecRed(pwfl_s2c, pwfl_s2c, p2r, true);
  // check the validity of s2c
  for(long i = 0; i < n_slots; i++) {
    NTL::ZZ val_in;
    if(NTL::deg(slots_in[i]) == -1)
      val_in = 0;
    else
      val_in = NTL::rem(slots_in[i][0], p2r);
    assertEq(val_in, pwfl_s2c[i * d], "something wrong with s2c?");
  }
  // try
  // for(long i = 0; i < n_slots; i++) {
  //   auto it = std::find(pwfl_s2c.begin(), pwfl_s2c.end(), NTL::ZZ(i+1));
  //   if(it != pwfl_s2c.end())
  //     printf("%ld found at %ld\n", i+1, it - pwfl_s2c.begin());
  // }
  // XXX: end debug
#endif

  // Mod-switch down if needed
  IndexSet s = ctxt.getPrimeSet() / context.getSpecialPrimes();
  assertTrue(s <= context.getCtxtPrimes(), "prime set is messed up");
  if (newKSFlag) {
    s.clear();
    s.insert(context.getIndexQks());
    ctxt.bringToSet(s);
  } else {
    if (s.card() > 3) { // leave only first three ciphertext primes
      long first = s.first();
      IndexSet s3(first, first + 2);
      s.retain(s3);
    }
    ctxt.modDownToSet(s);
  }


#ifdef HELIB_DEBUG
  CheckCtxt(ctxt, "after mod down to qKS");
  // XXX: debug
  NTL::ZZX poly_beforeKS;
  NTL::vec_ZZ pwfl_beforeKS;
  dbgKey->Decrypt(poly_beforeKS, ctxt);
  p2d_conv.ZZXtoPowerful(pwfl_beforeKS, poly_beforeKS);
  vecRed(pwfl_beforeKS, pwfl_beforeKS, p2r, true);
  for(long i = 0; i < n_slots; i++) {
    assertEq(pwfl_s2c[i * d], pwfl_beforeKS[i * d], "something wrong with mod down to qKS");
  }
  // XXX: end debug
#endif


  // key-switch to the bootstrapping key
  ctxt.reLinearize(recryptKeyID, newKSFlag);

#ifdef HELIB_DEBUG
  CheckCtxt(ctxt, "after key switching");
  // XXX: debug
  NTL::ZZX poly_beforeRaw;
  NTL::vec_ZZ pwfl_beforeRaw;
  dbgKey->Decrypt(poly_beforeRaw, ctxt);
  p2d_conv.ZZXtoPowerful(pwfl_beforeRaw, poly_beforeRaw);
  vecRed(pwfl_beforeRaw, pwfl_beforeRaw, p2r, true);
  for(long i = 0; i < n_slots; i++) {
    assertEq(pwfl_beforeKS[i * d], pwfl_beforeRaw[i * d], "something wrong with ks to btk");
  }
  // XXX: end debug
#endif


  // "raw mod-switch" to the bootstrapping mosulus q=p^e+1.
  std::vector<NTL::ZZX> zzParts; // the mod-switched parts, in ZZX format

  double mfac = ctxt.getContext().getZMStar().getNormBnd();
  double noise_est;
  if (newBtsFlag)
    noise_est = ctxt.rawModSwitchNew(zzParts, q) * mfac;
  else
    noise_est = ctxt.rawModSwitch(zzParts, q) * mfac;
#ifdef HELIB_DEBUG
  // XXX: debug, check the bound on I after native bts
  NTL::ZZX native_poly;
  NTL::vec_ZZ native_pwfl;
  rawDecrypt(native_poly, zzParts, dbgKey->getRecryptKey());
  p2d_conv.ZZXtoPowerful(native_pwfl, native_poly);
  auto native_pwfl_bound = NTL::conv<NTL::xdouble>(largestCoeff(native_pwfl)) /
                           NTL::conv<NTL::xdouble>(q);
  std::cout << "native pwfl bound is " << native_pwfl_bound << "\n";
  // XXX: end debug native
#endif

  // noise_est is an upper bound on the L-infty norm of the scaled noise
  // in the pwrfl basis
  double beta =
      p2r * (context.getScale() * sqrt(double(phim) / 12.0) *
                 (sqrt(double(context.getHwt()) * log(double(phim))) + 1) +
             0.5);
  double noise_bnd = newBtsFlag ? 2 * beta * mfac
                                : HELIB_MIN_CAP_FRAC * p2r *
                                      ctxt.getContext().boundForRecryption();
  // noise_bnd is the bound assumed in selecting the parameters
  double noise_rat = noise_est / noise_bnd;

  HELIB_STATS_UPDATE("raw-mod-switch-noise", noise_rat);

  if (noise_rat > 1) {
    // TODO: Turn the following preprocessor logics into a warnOrThrow function
    std::string message =
        "rawModSwitch scaled noise exceeds bound: " + std::to_string(noise_rat);
#ifdef HELIB_DEBUG
    Warning(message);
#else
    throw LogicError(message);
#endif
  }

  assertEq(zzParts.size(),
           (std::size_t)2,
           "Exactly 2 parts required for mod-switching in thin bootstrapping");

#ifdef HELIB_DEBUG
  if (dbgKey) {
    if (!newBtsFlag)
      checkRecryptBounds(zzParts,
                         dbgKey->getRecryptKey(),
                         ctxt.getContext(),
                         q);
  }
#endif

  std::vector<NTL::ZZX> v;
  v.resize(2);

  Ctxt modUpCtxt(ctxt.getPubKey(), p2eNew);
  std::vector<NTL::ZZX> I_part;

  // debug variables

  NTL::ZZX poly_before_mod;
  NTL::vec_ZZ pwfl_before_mod;
  NTL::ZZX poly_after_mod, I_poly, modUpPoly;
  NTL::vec_ZZ I_part_pwfl, pwfl_diff, pwfl_after_mod, modUpPwfl;
  std::vector<NTL::ZZX> slots_scaleUp;

  if (!newBtsFlag) {
    // Add multiples of q to make the zzParts divisible by p^{e'}
    for (long i : range(2)) {
      // make divisible by p^{e'}

      newMakeDivisible(zzParts[i], p2ePrime, q, ctxt.getContext(), v[i]);
    }
  } else {
#ifdef HELIB_DEBUG
    // check the capacity of zzParts
    rawDecrypt(poly_before_mod, zzParts, dbgKey->getRecryptKey());
    p2d_conv.ZZXtoPowerful(pwfl_before_mod, poly_before_mod);
    vecRed(pwfl_before_mod, pwfl_before_mod, q, false);
    NTL::ZZ maxPwfl = largestCoeff(pwfl_before_mod);
    double capacity =
        log(double(q) / NTL::to_long(maxPwfl) / 2.0 / aux) / log(2.0);
    printf("pwfl capacity of zzParts is %f\n", capacity);
    assertTrue(capacity > 0, "not enough capacity before scale");
    // check the correctness of raw MS
    long divFac = NTL::MulMod(q % p2r, intFactor, p2r);
    divFac = NTL::InvMod(divFac, p2r);
    for(long i = 0; i < n_slots; i++) {
      long before = NTL::to_long(pwfl_beforeRaw[i * d]);
      long after = NTL::rem(pwfl_before_mod[i * d], p2r);
      after = NTL::MulMod(after, divFac, p2r);
      assertEq(after, before, "something wrong with raw ms");
    }
#endif
    // NOTE: in thin bootstrapping, the Zpr values in slots are first
    //  move to the coeff domain by slotToCoeff
    //  then the multAndGetOverflowPart adds noise to the lower digits
    //  and the unoccupied coeffs
    //  a coeffToSlot map moves the noisy coeffs into the slots,
    //  then a trace-like map removes the values out of Z_p2eNew
    //  In other words, we have apply the second map to
    //  both modUpCtxt and Enc(I+aux*I', p2eNew)
    multAndGetOverflowPart(zzParts, I_part, ctxt.getContext());
    // NOTE: I_part is discarded

#ifdef HELIB_DEBUG
    // XXX: debug
    rawDecrypt(poly_after_mod, zzParts, dbgKey->getRecryptKey());
    // check the overflow part
    rawDecrypt(I_poly, I_part, dbgKey->getRecryptKey());
    p2d_conv.ZZXtoPowerful(I_part_pwfl, I_poly);
    vecRed(I_part_pwfl, I_part_pwfl, aux, false);
    p2d_conv.powerfulToZZX(I_poly, I_part_pwfl);

    p2d_conv.ZZXtoPowerful(pwfl_diff, poly_after_mod - I_poly * q);
    // vecRed(pwfl_diff, pwfl_diff, NTL::ZZ(aux) * NTL::ZZ(q), false);
    for (long i = 0; i < phim; i++) {
      assertTrue(bool(NTL::abs(pwfl_diff[i]) <= (q / 2)),
                 "overflow part wrong?");
    }
    // assertTrue(NTL::to_long(largestCoeff(pwfl_diff)) <= (q / 2),
    //            "overflow part wrong?");
    // check the relationship between pwfl_before_mod and pwfl_after_mod
    // now check if the bound on I is valid
    p2d_conv.ZZXtoPowerful(pwfl_after_mod, poly_after_mod);
    NTL::xdouble actualBoundI =
        NTL::conv<NTL::xdouble>(largestCoeff(pwfl_after_mod)) / NTL::xdouble(q);
    std::cerr << "actual bound on I is " << actualBoundI << "\n";
    assertTrue(
        bool(NTL::floor(actualBoundI) <= ceil(context.boundForRecryption())),
        "bound on I exceeded");
    vecRed(pwfl_after_mod, pwfl_after_mod, q, false);
    if (t < 0) {
      for (long i = 0; i < phim; i++) {
        long before = NTL::rem(pwfl_before_mod[i], p2r);
        before = NTL::MulMod(before, aux % p2r, p2r);
        long after = NTL::rem(pwfl_after_mod[i], p2r);
        assertEq(before,
                 after,
                 "something wrong with upscaling, non-power-of-p aux");
      }
    } else {
      for (long i = 0; i < phim; i++) {
        long before = NTL::rem(pwfl_before_mod[i], p2r);
        long after = NTL::rem(pwfl_after_mod[i], p2eNew);
        assertEq(before * aux,
                 after,
                 "something wrong with upscaling, power-of-p aux");
      }
    }
// XXX: end debug
#endif

    // now handle the mod-up ctxt
    modUpCtxt.ptxtSpace = p2eNew;
    modUpCtxt.primeSet = context.getCtxtPrimes() | context.getModUpPrimes();
    long QmodP = 1;
    // be careful with single-precision modular arithmetic...
    for (auto i : modUpCtxt.primeSet)
      QmodP = NTL::MulMod(QmodP, context.ithPrime(i) % p2eNew, p2eNew);
    // we pretend the IMPLIED INTFACTOR of m' mod P is 1
    // (i.e., zzParts after raw-ms decrypts to m')
    // which is actually intFactor * [q]_{p^r}
    // we want the mult-by-aux and modded-up ctxt 
    // to decrypt to p^t*m'+[q]_{p^eNew}*I...
    // aux = p^t, stored intFactor = [Q^-1]_P
    modUpCtxt.intFactor = NTL::InvMod(QmodP, p2eNew);
    // divide by [q]_p2eNew, so that the encrypted value is I+[q]^-1_{p2eNew}*p^t*m'
    modUpCtxt.intFactor = NTL::MulMod(modUpCtxt.intFactor, q % p2eNew, p2eNew);
    modUpCtxt.noiseBound = context.boundForRecryption() * q;
    for (int i = 0; i < 2; i++)
      modUpCtxt.addPart(DoubleCRT(zzParts[i], context, modUpCtxt.primeSet),
                        ctxt.parts[i].skHandle);
    // discard the mod-up primes
    modUpCtxt.modDownToSet(context.getCtxtPrimes());

    // switch to the normal sk (keyID = 0)
    modUpCtxt.reLinearize();

#ifdef HELIB_DEBUG
    dbgKey->Decrypt(modUpPoly, modUpCtxt);
    p2d_conv.ZZXtoPowerful(modUpPwfl, modUpPoly);
    vecRed(modUpPwfl, modUpPwfl, p2eNew, false);
    p2d_conv.ZZXtoPowerful(I_part_pwfl, I_poly);
    vecRed(I_part_pwfl, I_part_pwfl, aux, false);

    p2d_conv.ZZXtoPowerful(pwfl_after_mod, poly_after_mod);
    vecRed(pwfl_after_mod, pwfl_after_mod, q, false);

    // for aux = p^t
    // the message in zzParts befure "multAndGetOverflowPart" is m' mod p2r
    // after that, the scaled message becomes p^t*m' mod p2eNew
    // by modding up, the message becomes p^t*m' + [q]_p2eNew*I
    // setting the modUpCtxt.intFactor = [Q^-1]_p2eNew ensures
    // modUpCtxt decrypts to p^t*m' + [q]_p2eNew*I
    for (long i = 0; i < phim; i++) {
      if (t > 0) {
        assertEq(NTL::rem(modUpPwfl[i] - I_part_pwfl[i] * (q % p2eNew) -
                              pwfl_after_mod[i],
                          p2eNew),
                 0L,
                 "something wrong with assembly...? power-of-p aux");
      } else {
        // aux is not a power of p
        // old message: m' mod p2r
        // scaled message: aux*m' mod p2r (= pwfl_after_mod)
        // modding up: aux*m' + [q]_p2r*I
        // by setting stored intFactor to aux*[Q]_p2r:
        //  m' + [q*aux^-1]_p2r*I
        assertEq(NTL::rem(aux * modUpPwfl[i] - I_part_pwfl[i] * (q % p2r) -
                              pwfl_after_mod[i],
                          p2r),
                 0L,
                 "something wrong with assembly...? non-power-of-p aux");
      }
    }
#endif
//     // NOTE: this map is exclusive to thin-bts
//     time_linear2_start_second = steady_clock::now();
//     trcData.coeffToSlot->apply(modUpCtxt);
//     time_linear2_end_second = steady_clock::now();
//     // std::cout << "time for 2nd linear2 is " << 
//     //   duration_cast<duration<double>>(time_linear2_end_second - time_linear2_start_second).count()
//     //   << "\n";
// #ifdef HELIB_DEBUG
//     ea_boot->decrypt(modUpCtxt, *dbgKey, slots_scaleUp);
// #endif
  }

#ifdef HELIB_DEBUG
  if (dbgKey) {
    if (!newBtsFlag) {
      checkRecryptBounds_v(v, dbgKey->getRecryptKey(), ctxt.getContext(), q);
      checkCriticalValue(zzParts,
                         dbgKey->getRecryptKey(),
                         ctxt.getContext().getRcData(),
                         q);
    }
  }
#endif

  if (!newBtsFlag)
    for (long i : range(zzParts.size())) {
      zzParts[i] /= p2ePrime; // divide by p^{e'}
    }

  // NOTE: here we lose the intFactor associated with ctxt.
  // We will restore it below.
  ctxt = recryptEkey;

  NTL::ZZX I_poly_new;
  NTL::vec_ZZ I_poly_new_pwfl;
  if (!newBtsFlag) {
    ctxt.multByConstant(zzParts[1]);
    ctxt.addConstant(zzParts[0]);
    benchmarker.bits_after_inner_prod = ctxt.capacity();
  } else {
    ctxt = modUpCtxt;
    benchmarker.bits_after_inner_prod = ctxt.capacity();
  }
#ifdef HELIB_DEBUG
  CheckCtxt(ctxt, "after bootKeySwitch");
#endif

  HELIB_NTIMER_STOP(AAA_bootKeySwitch);

  // Move the powerful-basis coefficients to the plaintext slots
  HELIB_NTIMER_START(AAA_coeffToSlot);
  auto bits_down_linear2 = ctxt.capacity();
  auto time_linear2_start = steady_clock::now();
  trcData.coeffToSlot->apply(ctxt);
  auto time_linear2_end = steady_clock::now();
  bits_down_linear2 -= ctxt.capacity();
  HELIB_NTIMER_STOP(AAA_coeffToSlot);
  // std::cout << "time for linear2 is " << 

#ifdef HELIB_DEBUG
  std::vector<NTL::ZZX> slots_I;
  ea_boot->decrypt(ctxt, *dbgKey, slots_I);
  // XXX: debug, check if the slots match
  auto vals_in = balrem_vec_ZZX(slots_in, p2r);
  auto vals_scaleUp = balrem_vec_ZZX(slots_scaleUp, p2eNew);
  auto val_I = balrem_vec_ZZX(slots_I, p2eNew);
  if (t > 0)
    for (long i = 0; i < n_slots; i++) {
      long tmp = balRem(NTL::rem(val_I[i], aux), aux); // I without I*
      tmp = (tmp + p2eNew) % p2eNew;
      tmp = NTL::MulMod(q % p2eNew, tmp, p2eNew);
      tmp = NTL::SubMod(NTL::rem(vals_scaleUp[i], p2eNew), tmp, p2eNew);
      assertEq(tmp % aux, 0L, "not divisible by aux?");
      tmp /= aux;
      tmp %= p2r;
      long factorDiv = intFactor;
      factorDiv = NTL::MulMod(factorDiv % p2r, q % p2r, p2r);
      assertEq(tmp,
               NTL::MulMod(NTL::rem(vals_in[i], p2r), factorDiv, p2r),
               "not matching, p^t");
    }
    // TODO: t < 0
    // XXX: end debug
#endif

#ifdef HELIB_DEBUG
  CheckCtxt(ctxt, "after coeffToSlot");
#endif

  // Extract the digits e-e'+r-1,...,e-e' (from fully packed slots)
  HELIB_NTIMER_START(AAA_extractDigitsThin);
  auto bits_down_extract = ctxt.capacity();
  auto time_extract_start = steady_clock::now();
  extractDigitsThin(ctxt, e - ePrime, r, ePrime, true);
  auto time_extract_end = steady_clock::now();
  bits_down_extract -= ctxt.capacity();
  HELIB_NTIMER_STOP(AAA_extractDigitsThin);

#ifdef HELIB_DEBUG
  if (newBtsFlag) {
    std::vector<NTL::ZZX> slots_I_out;
    ea_boot->decrypt(ctxt, *dbgKey, slots_I_out);
    auto vals_I_out = balrem_vec_ZZX(slots_I_out, p2eNew);
    for (long i = 0; i < n_slots; i++) {
      assertEq(balRem(NTL::rem(val_I[i], aux), aux),
               NTL::to_long(vals_I_out[i]),
               "something wrong with digit extraction");
    }
  }
#endif

#ifdef HELIB_DEBUG
  CheckCtxt(ctxt, "after extractDigitsThin");
#endif

  // restore intFactor
  if (newBtsFlag) {
    // the ctxt encrypts [q]^-1_p2r*m' mod p2r
    ctxt *= q % p2r;
    // now correct the intFactor of m'
    // NOTE: for old bts, since q=p^e+1, [q]_p^r==1
    ctxt.intFactor =
        NTL::MulMod(ctxt.intFactor,
                    NTL::MulMod(intFactor, q % ptxtSpace, ptxtSpace),
                    ptxtSpace);
  } else if (intFactor != 1)
    ctxt.intFactor = NTL::MulMod(ctxt.intFactor, intFactor, ptxtSpace);
  auto time_boot_end = steady_clock::now();
  benchmarker.bits_down_extract = bits_down_extract;
  benchmarker.bits_down_linear_1 = bits_down_linear1;
  benchmarker.bits_down_linear_2 = bits_down_linear2;
  benchmarker.bits_final = ctxt.capacity();
  benchmarker.time_linear_1 =
      duration_cast<duration<double>>(time_linear1_end - time_linear1_start)
          .count();
  benchmarker.time_linear_2 =
      duration_cast<duration<double>>(time_linear2_end - time_linear2_start)
          .count();
  benchmarker.time_extract =
      duration_cast<duration<double>>(time_extract_end - time_extract_start)
          .count();
  benchmarker.time_total =
      duration_cast<duration<double>>(time_boot_end - time_boot_start).count();
  return benchmarker;
}


#ifdef HELIB_DEBUG

static void checkCriticalValue(const std::vector<NTL::ZZX>& zzParts,
                               const DoubleCRT& sKey,
                               const RecryptData& rcData,
                               long q)
{
  NTL::ZZX ptxt;
  rawDecrypt(ptxt, zzParts, sKey); // no mod q

  NTL::Vec<NTL::ZZ> powerful;
  rcData.p2dConv->ZZXtoPowerful(powerful, ptxt);
  NTL::xdouble max_pwrfl = NTL::conv<NTL::xdouble>(largestCoeff(powerful));
  double critical_value = NTL::conv<double>((max_pwrfl / q) / q);

  vecRed(powerful, powerful, q, false);
  max_pwrfl = NTL::conv<NTL::xdouble>(largestCoeff(powerful));
  critical_value += NTL::conv<double>(max_pwrfl / q);

  HELIB_STATS_UPDATE("critical-value", critical_value);

  std::cerr << "=== critical_value=" << critical_value;
  if (critical_value > 0.5)
    std::cerr << " BAD-BOUND";

  std::cerr << "\n";
}

static void checkRecryptBounds(const std::vector<NTL::ZZX>& zzParts,
                               const DoubleCRT& sKey,
                               const Context& context,
                               long q)
{
  const RecryptData& rcData = context.getRcData();
  double coeff_bound = context.boundForRecryption();
  long p2r = context.getAlMod().getPPowR();

  NTL::ZZX ptxt;
  rawDecrypt(ptxt, zzParts, sKey); // no mod q

  NTL::Vec<NTL::ZZ> powerful;
  rcData.p2dConv->ZZXtoPowerful(powerful, ptxt);
  double max_pwrfl = NTL::conv<double>(largestCoeff(powerful));
  double ratio = max_pwrfl / (2 * q * coeff_bound);

  HELIB_STATS_UPDATE("|x|/bound", ratio);

  std::cerr << "=== |x|/bound=" << ratio;
  if (ratio > 1.0)
    std::cerr << " BAD-BOUND";

  vecRed(powerful, powerful, q, false);
  max_pwrfl = NTL::conv<double>(largestCoeff(powerful));
  ratio = max_pwrfl / (2 * p2r * coeff_bound);

  HELIB_STATS_UPDATE("|x%q|/bound", ratio);

  std::cerr << ", (|x%q|)/bound=" << ratio;
  if (ratio > 1.0)
    std::cerr << " BAD-BOUND";

  std::cerr << "\n";
}

static void checkRecryptBounds_v(const std::vector<NTL::ZZX>& v,
                                 const DoubleCRT& sKey,
                                 const Context& context,
                                 UNUSED long q)
{
  const RecryptData& rcData = context.getRcData();

  long p = context.getP();
  long e = rcData.e;
  long p2e = NTL::power_long(p, e);
  long ePrime = rcData.ePrime;
  long p2ePrime = NTL::power_long(p, ePrime);
  long phim = context.getPhiM();

  double fudge = compute_fudge(p2ePrime, p2e);

  double coeff_bound = context.boundForRecryption() * fudge;

  double sigma = context.stdDevForRecryption() * fudge;

  NTL::ZZX ptxt;
  rawDecrypt(ptxt, v, sKey); // no mod q

  NTL::Vec<NTL::ZZ> powerful;
  rcData.p2dConv->ZZXtoPowerful(powerful, ptxt);
  double max_pwrfl = NTL::conv<double>(largestCoeff(powerful));

  double denom = p2ePrime * coeff_bound;
  double ratio = max_pwrfl / denom;

  HELIB_STATS_UPDATE("|v|/bound", ratio);

  std::cerr << "=== |v|/bound=" << ratio;
  if (ratio > 1.0)
    std::cerr << " BAD-BOUND";
  std::cerr << "\n";

  ptxt -= v[0]; // so now ptxt is just sKey * v[1]
  rcData.p2dConv->ZZXtoPowerful(powerful, ptxt);

  assertEq(powerful.length(), phim, "length should be phim");

  double ran_pwrfl = NTL::conv<double>(powerful[NTL::RandomBnd(phim)]);
  // pick a random coefficient in the poweful basis

  double std_devs = fabs(ran_pwrfl) / (p2ePrime * sigma);
  // number of standard deviations away from mean

  // update various indicator variables
  HELIB_STATS_UPDATE("sigma_0_5", double(std_devs <= 0.5)); // 0.383
  HELIB_STATS_UPDATE("sigma_1_0", double(std_devs <= 1.0)); // 0.683
  HELIB_STATS_UPDATE("sigma_1_5", double(std_devs <= 1.5)); // 0.866
  HELIB_STATS_UPDATE("sigma_2_0", double(std_devs <= 2.0)); // 0.954
  HELIB_STATS_UPDATE("sigma_2_5", double(std_devs <= 2.5)); // 0.988
  HELIB_STATS_UPDATE("sigma_3_0", double(std_devs <= 3.0)); // 0.997, 1 in 370
  HELIB_STATS_UPDATE("sigma_3_5",
                     double(std_devs <= 3.5)); // 0.999535, 1 in 2149
  HELIB_STATS_UPDATE("sigma_4_0",
                     double(std_devs <= 4.0)); // 0.999937, 1 in 15787

  // compute sample variance, and scale by the variance we expect
  HELIB_STATS_UPDATE("sigma_calc",
                     fsquare(ran_pwrfl) / fsquare(p2ePrime * sigma));

  // save the scaled value for application of other tests
  HELIB_STATS_SAVE("v_values", ran_pwrfl / (p2ePrime * sigma));
}

#endif

#if 0
void fhe_stats_print(long iter, const Context& context)
{
   long phim = context.getPhiM();

   std::cerr << "||||| recryption stats ||||\n";
   std::cerr << "**** averages ****\n";
   std::cerr << "=== critical_value=" << (fhe_stats_cv_sum/iter) << "\n";
   std::cerr << "=== |x|/bound=" << (fhe_stats_x_sum/iter) << "\n";
   std::cerr << "=== |x%q|/bound=" << (fhe_stats_xmod_sum/iter) << "\n";
   std::cerr << "=== |u|/bound=" << (fhe_stats_u_sum/iter) << "\n";
   std::cerr << "=== |v|/bound=" << (fhe_stats_v_sum/iter) << "\n";
   std::cerr << "**** maxima ****\n";
   std::cerr << "=== critical_value=" << (fhe_stats_cv_max) << "\n";
   std::cerr << "=== |x|/bound=" << (fhe_stats_x_max) << "\n";
   std::cerr << "=== |x%q|/bound=" << (fhe_stats_xmod_max) << "\n";
   std::cerr << "=== |u|/bound=" << (fhe_stats_u_max) << "\n";
   std::cerr << "=== |v|/bound=" << (fhe_stats_v_max) << "\n";
   std::cerr << "**** theoretical bounds ***\n";
   std::cerr << "=== single-max=" << (sqrt(2.0*log(phim))/context.scale) << "\n";
   std::cerr << "=== global-max=" << (sqrt(2.0*(log(iter)+log(phim)))/context.scale) << "\n";


}
#endif

} // namespace helib
