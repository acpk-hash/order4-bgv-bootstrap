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
#include <helib/Context.h>
#include <helib/polyEval.h>
#include <map>
#include <cmath>

namespace helib {

// Returns the e'th power of X, computing it as needed
Ctxt& DynamicCtxtPowers::getPower(long e)
{
  if (v.at(e - 1).isEmpty()) { // Not computed yet, compute it now

    // largest power of two smaller than e
    long k = 1L << (NTL::NextPowerOfTwo(e) - 1);
    v[e - 1] = getPower(e - k); // compute X^e = X^{e-k} * X^k
    v[e - 1].multiplyBy(getPower(k));
    // FIXME: could drop down / cleanup further as an optimization?
  }
  return v[e - 1];
}

// Local functions for polynomial evaluation in some special cases
static void simplePolyEval(Ctxt& ret,
                           const NTL::ZZX& poly,
                           DynamicCtxtPowers& babyStep);
static void PatersonStockmeyer(Ctxt& ret,
                               const NTL::ZZX& poly,
                               long k,
                               long t,
                               long delta,
                               DynamicCtxtPowers& babyStep,
                               DynamicCtxtPowers& giantStep);
static void degPowerOfTwo(Ctxt& ret,
                          const NTL::ZZX& poly,
                          long k,
                          DynamicCtxtPowers& babyStep,
                          DynamicCtxtPowers& giantStep);
static void recursivePolyEval(Ctxt& ret,
                              const NTL::ZZX& poly,
                              long k,
                              DynamicCtxtPowers& babyStep,
                              DynamicCtxtPowers& giantStep);

static void recursivePolyEval(Ctxt& ret,
                              const Ctxt poly[],
                              long nCoeffs,
                              const NTL::Vec<Ctxt>& powers);

static void recursivePolyEvalNew(Ctxt& ret,
                                 const NTL::ZZX& poly,
                                 long k,
                                 DynamicCtxtPowers& babyStep,
                                 std::vector<Ctxt>& giantStep);

// Main entry point: Evaluate an encrypted polynomial on an encrypted input
// return in ret = sum_i poly[i] * x^i
void polyEval(Ctxt& ret, const NTL::Vec<Ctxt>& poly, const Ctxt& x)
{
  if (poly.length() <= 1) { // Some special cases
    if (poly.length() == 0)
      ret.clear(); // empty polynomial
    else
      ret = poly[0]; // constant polynomial
    return;
  }
  long deg = poly.length() - 1;

  long logD = NTL::NextPowerOfTwo(divc(poly.length(), 3));
  long d = 1L << logD;

  // We have d <= deg(poly) < 3d
  assertInRange(deg, d, 3l * d, "Poly degree not in [d, 3d)");

  NTL::Vec<Ctxt> powers(NTL::INIT_SIZE, logD + 1, x);
  if (logD > 0) {
    powers[1].square();
    for (long i = 2; i <= logD; i++) { // powers[i] = x^{2^i}
      powers[i] = powers[i - 1];
      powers[i].square();
    }
  }

  // Compute in three parts p0(X) + ( p1(X) + p2(X)*X^d )*X^d
  Ctxt tmp(ZeroCtxtLike, ret);
  recursivePolyEval(ret,
                    &poly[d],
                    std::min(d, poly.length() - d),
                    powers); // p1(X)

  if (poly.length() > 2 * d) { // p2 is not empty
    recursivePolyEval(tmp,
                      &poly[2 * d],
                      poly.length() - 2 * d,
                      powers); // p2(X)
    tmp.multiplyBy(powers[logD]);
    ret += tmp;
  }
  ret.multiplyBy(powers[logD]); // ( p1(X) + p2(X)*X^d )*X^d

  recursivePolyEval(tmp, &poly[0], d, powers); // p0(X)
  ret += tmp;
}

static void recursivePolyEval(Ctxt& ret,
                              const Ctxt poly[],
                              long nCoeffs,
                              const NTL::Vec<Ctxt>& powers)
{
  if (nCoeffs <= 1) { // edge condition
    if (nCoeffs == 0)
      ret.clear(); // empty polynomial
    else
      ret = poly[0]; // constant polynomial
    return;
  }
  long logD = NTL::NextPowerOfTwo(nCoeffs) - 1;
  long d = 1L << logD;
  Ctxt tmp(ZeroCtxtLike, ret);
  recursivePolyEval(tmp, &(poly[d]), nCoeffs - d, powers);
  recursivePolyEval(ret, &(poly[0]), d, powers);
  tmp.multiplyBy(powers[logD]);
  ret += tmp;
}

// Main entry point: Evaluate a cleartext polynomial on an encrypted input
void polyEval(Ctxt& ret, NTL::ZZX poly, const Ctxt& x, long k)
// Note: poly is passed by value, so caller keeps the original
{
  if (deg(poly) <= 2) {  // nothing to optimize here
    if (deg(poly) < 1) { // A constant
      ret.clear();
      ret.addConstant(coeff(poly, 0));
    } else { // A linear or quadratic polynomial
      DynamicCtxtPowers babyStep(x, deg(poly));
      simplePolyEval(ret, poly, babyStep);
    }
    return;
  }

  // How many baby steps: set k~sqrt(n/2), rounded up/down to a power of two

  // FIXME: There may be some room for optimization here: it may be possible
  // to choose k as something other than a power of two and still maintain
  // optimal depth, in principle we can try all possible values of k between
  // two consecutive powers of two and choose the one that gives the least
  // number of multiplies, conditioned on minimum depth.

  if (k <= 0) {
    long kk = (long)sqrt(deg(poly) / 2.0);
    k = 1L << NTL::NextPowerOfTwo(kk);

    // heuristic: if k>>kk then use a smaller power of two
    if ((k == 16 && deg(poly) > 167) || (k > 16 && k > (1.44 * kk)))
      k /= 2;
  }
#ifdef HELIB_DEBUG
  std::cerr << "  k=" << k;
#endif

  long n = divc(deg(poly), k); // n = ceil(deg(p)/k), deg(p) >= k*n
  DynamicCtxtPowers babyStep(x, k);
  const Ctxt& x2k = babyStep.getPower(k);

  // Special case when deg(p)>k*(2^e -1)
  if (n == (1L << NTL::NextPowerOfTwo(n))) { // n is a power of two
    DynamicCtxtPowers giantStep(x2k, n / 2);
    degPowerOfTwo(ret, poly, k, babyStep, giantStep);
    return;
  }

  // If n is not a power of two, ensure that poly is monic and that
  // its degree is divisible by k, then call the recursive procedure

  const NTL::ZZ p = NTL::to_ZZ(x.getPtxtSpace());
  NTL::ZZ top = LeadCoeff(poly);
  NTL::ZZ topInv; // the inverse mod p of the top coefficient of poly (if any)
  bool divisible = (n * k == deg(poly)); // is the degree divisible by k?
  long nonInvertible = InvModStatus(topInv, top, p);
  // 0 if invertible, 1 if not

  // FIXME: There may be some room for optimization below: instead of
  // adding a term X^{n*k} we can add X^{n'*k} for some n'>n, so long
  // as n' is smaller than the next power of two. We could save a few
  // multiplications since giantStep[n'] may be easier to compute than
  // giantStep[n] when n' has fewer 1's than n in its binary expansion.

  // extra!=0 denotes an added term extra*X^{n*k}
  NTL::ZZ extra = NTL::ZZ::zero();
  if (!divisible || nonInvertible) { // need to add a term
    top = NTL::to_ZZ(1);             // new top coefficient is one
    topInv = top;                    // also the new inverse is one
    // set extra = 1 - current-coeff-of-X^{n*k}
    extra = SubMod(top, coeff(poly, n * k), p);
    SetCoeff(poly, n * k); // set the top coefficient of X^{n*k} to one
  }

  long t = IsZero(extra) ? divc(n, 2) : n;
  DynamicCtxtPowers giantStep(x2k, t);

  if (!IsOne(top)) {
    poly *= topInv; // Multiply by topInv to make into a monic polynomial
    for (long i = 0; i <= n * k; i++)
      rem(poly[i], poly[i], p);
    poly.normalize();
  }
  recursivePolyEval(ret, poly, k, babyStep, giantStep);

  if (!IsOne(top)) {
    ret.multByConstant(top);
  }

  if (!IsZero(extra)) { // if we added a term, now is the time to subtract back
    Ctxt topTerm = giantStep.getPower(n);
    topTerm.multByConstant(extra);
    ret -= topTerm;
  }
}

// inputs are always relined, outputs are not
long cost_n_leaves(long n) {
  long cost = 0;
  if(n <= 2)
    return 0;  
  long prevPo2 = NTL::NextPowerOfTwo(n);
  if(1L << prevPo2 > n)
    prevPo2 -= 1;
  prevPo2 = 1 << prevPo2;
  // a subtree of size 2^k has a cost of 2^(k-1)-2
  if(prevPo2 == n)
    return prevPo2 / 2 - 1;
  // 1 = relin the one the inputs
  cost += 1 + cost_n_leaves(n - prevPo2) + prevPo2 / 2 - 1;
  return cost;
}

void polyEvalNew(std::vector<Ctxt*>& ret,
                 NTL::vec_ZZX polys,
                 const Ctxt& x,
                 long k)
{
  std::vector<long> polyDegs(polys.length());
  long max_degree = 0;
  // make these polys odd
  for (long _ = 0; _ < polys.length(); _++) {
    auto& poly = polys[_];
    // note: moved to buildBtsPolys in extractDigits.cpp
    // for (long i = 0; i <= deg(poly); i+=2) {
    //   NTL::SetCoeff(poly, i, 0L);
    // }
    // poly.normalize();
    polyDegs[_] = deg(poly);
    max_degree = std::max(max_degree, deg(poly));
  }
  // compute the baby step size k, defaults to sqrt(max_degree)
  k = round(sqrt(max_degree + 1));
  // find the k with lowest cost
  long hi = 1L << NTL::NextPowerOfTwo(k);
  long lo = std::max(2L, hi >> 2);
  std::vector<long> minDepths(polyDegs.size());
  for (size_t i = 0; i < polyDegs.size(); i++)
    minDepths[i] = ceil(log(polyDegs[i]) / log(2));
  long minCost = max_degree * polyDegs.size();
  for (long kk = lo; kk <= hi; kk += 2) { // only search for even k
    bool depthOk = true;
    // shared computation, building the basis
    long cost = kk / 2 + ceil(log(max_degree + 1) / log(2.0));
    for (size_t j = 0; j < polyDegs.size(); j++) {
      long l = ceil(log(double(polyDegs[j] + 1) / kk) / log(2));
      long depth = l + ceil(log(kk) / log(2));
      if (minDepths[j] < depth) {
        depthOk = false;
        break; // we want optimal depth
      }
      long n_leaves = (polyDegs[j] + 1 + kk - 1) / kk;
      cost += cost_n_leaves(n_leaves);
    }
    if (!depthOk)
      continue;
    if (cost < minCost) {
      minCost = cost;
      k = kk;
    }
  }
#ifdef HELIB_DEBUG
  std::cerr << "  k=" << k;
#endif

  // now start the bsgs computation
  // first precompute the baby-step and giant-step polynomials
  long l = ceil(log(double(max_degree + 1) / k) / log(2));
  DynamicCtxtPowers babyStep(x, k);
  const Ctxt& x2k = babyStep.getPower(k);
  std::vector<Ctxt> x2kPow2(
      l,
      Ctxt(x.getPubKey(),
           x.getPtxtSpace())); // x^(2^i*k), for i = 0, 1, ..., l-1
  x2kPow2[0] = x2k;
  for (long i = 1; i < l; i++) {
    x2kPow2[i] = x2kPow2[i - 1];
    x2kPow2[i].square();
  }
  // now recurse
  for (long i = 0; i < polys.length(); i++) {
    recursivePolyEvalNew(*ret[i], polys[i], k, babyStep, x2kPow2);
    ret[i]->reLinearize();
  }
}

// Simple evaluation sum f_i * X^i, assuming that babyStep has enough powers
static void simplePolyEval(Ctxt& ret,
                           const NTL::ZZX& poly,
                           DynamicCtxtPowers& babyStep)
{
  ret.clear();
  if (deg(poly) < 0)
    return; // the zero polynomial always returns zero

  // Ensure that we have enough powers

  // ensure that we have enough powers
  assertTrue(deg(poly) <= babyStep.size(),
             "BabyStep has not enough powers "
             "(required more than deg(poly))");

  NTL::ZZ coef;
  NTL::ZZ p = NTL::to_ZZ(babyStep[0].getPtxtSpace());
  for (long i = 1; i <= deg(poly); i++) {
    rem(coef, coeff(poly, i), p);
    if (coef > p / 2)
      coef -= p;

    Ctxt tmp = babyStep.getPower(i); // X^i
    tmp.multByConstant(coef);        // f_i X^i
    ret += tmp;
  }
  // Add the free term
  rem(coef, ConstTerm(poly), p);
  if (coef > p / 2)
    coef -= p;
  ret.addConstant(coef);
  //  if (verbose) checkPolyEval(ret, babyStep[0], poly);
}

// The recursive procedure in the Paterson-Stockmeyer
// polynomial-evaluation algorithm from SIAM J. on Computing, 1973.
// This procedure assumes that poly is monic, deg(poly)=k*(2t-1)+delta
// with t=2^e, and that babyStep contains >= k+delta powers
static void PatersonStockmeyer(Ctxt& ret,
                               const NTL::ZZX& poly,
                               long k,
                               long t,
                               long delta,
                               DynamicCtxtPowers& babyStep,
                               DynamicCtxtPowers& giantStep)
{
  if (deg(poly) <= babyStep.size()) { // Edge condition, use simple eval
    simplePolyEval(ret, poly, babyStep);
    return;
  }
  NTL::ZZX r = trunc(poly, k * t);      // degree <= k*2^e-1
  NTL::ZZX q = RightShift(poly, k * t); // degree == k(2^e-1) +delta

  const NTL::ZZ p = NTL::to_ZZ(babyStep[0].getPtxtSpace());
  const NTL::ZZ& coef = coeff(r, deg(q));
  SetCoeff(r, deg(q), coef - 1); // r' = r - X^{deg(q)}

  NTL::ZZX c, s;
  DivRem(c, s, r, q); // r' = c*q + s
  // deg(s)<deg(q), and if c!= 0 then deg(c)<k-delta

  assertTrue(deg(s) < deg(q), "Degree of s is not less than degree of q");
  assertTrue(IsZero(c) || deg(c) < k - delta,
             "Nonzero c has not degree smaller than k - delta");
  SetCoeff(s, deg(q)); // s' = s + X^{deg(q)}, deg(s)==deg(q)

  // reduce the coefficients modulo p
  for (long i = 0; i <= deg(c); i++)
    rem(c[i], c[i], p);
  c.normalize();
  for (long i = 0; i <= deg(s); i++)
    rem(s[i], s[i], p);
  s.normalize();

  // Evaluate recursively poly = (c+X^{kt})*q + s'
  PatersonStockmeyer(ret, q, k, t / 2, delta, babyStep, giantStep);

  Ctxt tmp(ret.getPubKey(), ret.getPtxtSpace());
  simplePolyEval(tmp, c, babyStep);
  tmp += giantStep.getPower(t);
  ret.multiplyBy(tmp);

  PatersonStockmeyer(tmp, s, k, t / 2, delta, babyStep, giantStep);
  ret += tmp;
}

// This procedure assumes that k*(2^e +1) > deg(poly) > k*(2^e -1),
// and that babyStep contains >= k + (deg(poly) mod k) powers
static void degPowerOfTwo(Ctxt& ret,
                          const NTL::ZZX& poly,
                          long k,
                          DynamicCtxtPowers& babyStep,
                          DynamicCtxtPowers& giantStep)
{
  if (deg(poly) <= babyStep.size()) { // Edge condition, use simple eval
    simplePolyEval(ret, poly, babyStep);
    return;
  }
  long n = deg(poly) / k;                     // We assume n=2^e or n=2^e -1
  n = 1L << NTL::NextPowerOfTwo(n);           // round up to n=2^e
  NTL::ZZX r = trunc(poly, (n - 1) * k);      // degree <= k(2^e-1)-1
  NTL::ZZX q = RightShift(poly, (n - 1) * k); // 0 < degree < 2k
  SetCoeff(r, (n - 1) * k);                   // monic, degree == k(2^e-1)
  q -= 1;

  PatersonStockmeyer(ret, r, k, n / 2, 0, babyStep, giantStep);

  Ctxt tmp(ret.getPubKey(), ret.getPtxtSpace());
  simplePolyEval(tmp, q, babyStep); // evaluate q

  // multiply by X^{k(n-1)} with minimum depth
  for (long i = 1; i < n; i *= 2) {
    tmp.multiplyBy(giantStep.getPower(i));
  }
  ret += tmp;
}

static void recursivePolyEval(Ctxt& ret,
                              const NTL::ZZX& poly,
                              long k,
                              DynamicCtxtPowers& babyStep,
                              DynamicCtxtPowers& giantStep)
{
  if (deg(poly) <= babyStep.size()) { // Edge condition, use simple eval
    simplePolyEval(ret, poly, babyStep);
    return;
  }

  long delta = deg(poly) % k;              // deg(poly) mod k
  long n = divc(deg(poly), k);             // ceil( deg(poly)/k )
  long t = 1L << (NTL::NextPowerOfTwo(n)); // t >= n, so t*k >= deg(poly)

  // Special case for deg(poly) = k * 2^e +delta
  if (n == t) {
    degPowerOfTwo(ret, poly, k, babyStep, giantStep);
    return;
  }

  // When deg(poly) = k*(2^e -1) we use the Paterson-Stockmeyer recursion
  if (n == t - 1 && delta == 0) {
    PatersonStockmeyer(ret, poly, k, t / 2, delta, babyStep, giantStep);
    return;
  }

  t = t / 2;

  // In any other case we have kt < deg(poly) < k(2t-1). We then set
  // u = deg(poly) - k*(t-1) and poly = q*X^u + r with deg(r)<u
  // and recurse on poly = (q-1)*X^u + (X^u+r)

  long u = deg(poly) - k * (t - 1);
  NTL::ZZX r = trunc(poly, u);      // degree <= u-1
  NTL::ZZX q = RightShift(poly, u); // degree == k*(t-1)
  q -= 1;
  SetCoeff(r, u); // degree == u

  PatersonStockmeyer(ret, q, k, t / 2, 0, babyStep, giantStep);

  Ctxt tmp = giantStep.getPower(u / k);
  if (delta != 0) { // if u is not divisible by k then compute it
    tmp.multiplyBy(babyStep.getPower(delta));
  }
  ret.multiplyBy(tmp);

  recursivePolyEval(tmp, r, k, babyStep, giantStep);
  ret += tmp;
}

static void recursivePolyEvalNew(Ctxt& ret,
                                 const NTL::ZZX& poly,
                                 long k,
                                 DynamicCtxtPowers& babyStep,
                                 std::vector<Ctxt>& giantStep)
{
  if (deg(poly) <= babyStep.size()) {
    simplePolyEval(ret, poly, babyStep);
    return;
  }
  long polyDeg = deg(poly);
  // find minimum i s.t. 2^(i+1)*k >= polyDeg + 1, then degSplit = 2^i*k
  // i = ceil(log_2((polyDeg + 1) / k)) - 1
  long i = long(ceil(log(double(polyDeg + 1) / k) / log(2))) - 1;
  long degSplit = k << i;
  NTL::ZZX lowerHalf = trunc(poly, degSplit);
  NTL::ZZX higherHalf = RightShift(poly, degSplit);
  Ctxt tmp(babyStep.getPower(1).getPubKey(),
           babyStep.getPower(1).getPtxtSpace());
  recursivePolyEvalNew(tmp, lowerHalf, k, babyStep, giantStep);
  recursivePolyEvalNew(ret, higherHalf, k, babyStep, giantStep);
  // only relinearize before mult
  // the total number of relin is 2^(l-1)
  ret.reLinearize();
  ret.multLowLvl(giantStep[i]);
  ret += tmp;
}

// raise ciphertext to some power
void Ctxt::power(long e)
{
  if (e < 1)
    throw InvalidArgument("Cannot raise a ctxt to a non positive exponent");

  if (e == 1)
    return; // nothing to do

  long ell = NTL::NumBits(e); // e < 2^l <= 2e

  if (static_cast<unsigned long>(e) ==
      (1UL << (ell - 1))) { // e is a power of two, just square enough times
    while (--ell > 0)
      square();
    return;
  }

  // Otherwise use the "DynamicCtxtPowers" from polyEval, it uses e Ctxt
  // objects as temporary space but keeps the level as low as possible
  DynamicCtxtPowers pwrs(*this, e);
  *this = pwrs.getPower(e);
}

#if 0
/**********************************************************************/
/*     FOR DEBUGGING PURPOSES, the same procedure for plaintext x     */
/**********************************************************************/

static long verbose = 0;
class DynamicPtxtPowers {
private:
  long p;  // the modulus
  vec_long v;    // A vector storing the powers themselves
  vec_long dpth; // keeping track of depth for debugging purposes

public:
  DynamicPtxtPowers(long _x, long _p, long nPowers, long _d=1) : p(_p)
  {
    // Sanity check
    assertTrue<InvalidArgument>(_x >= 0l, "_x must be greater equal than 0");
    assertTrue<InvalidArgument>(_p > 1l, "_p must be greater than 1"); // Sanity check
    assertTrue<InvalidArgument>(nPowers > 0l, "nPowers must be greater than 0"); // Sanity check
    v.SetLength(nPowers);
    dpth.SetLength(nPowers);
    for (long i=1; i<nPowers; i++) // Initializes nPowers empty slots
      dpth[i]=v[i]=-1;

    v[0] = _x % p;         // store X itself in v[0]
    dpth[0] = _d;
    //    std::cerr << " initial power="<<v[0]<<", initial depth="<<dpth[0];
  }

  // Returns the e'th power, computing it as needed
  long getPower(long e); // must use e >= 1, else throws an exception

  // dp.at(i) and dp[i] both return the i+1st power
  long at(long i) { return getPower(i+1); }
  long operator[](long i) { return getPower(i+1); }

  const vec_long& getVector() const { return v; }
  long size() const { return v.length(); }
  bool wasComputed(long i) { return (i>=0 && i<v.length() && v[i]>=0); }
  long getDepth(long e) { return (e>0)? dpth.at(e-1) : 0; }
};

// Returns the e'th power, computing it as needed
long DynamicPtxtPowers::getPower(long e)
{
  // FIXME: Do we want to allow the vector to grow? If so then begin by
  // checking e<v.length() and resizing if not. Currently throws an exception.

  if (v.at(e-1)<0) { // Not computed yet, compute it now

    long k = 1L<<(NTL::NextPowerOfTwo(e)-1); // largest power of two smaller than e

    v[e-1] = getPower(e-k);             // compute X^e = X^{e-k} * X^k
    v[e-1] = MulMod(v[e-1], getPower(k), p);
    dpth[e-1] = max(getDepth(k),getDepth(e-k)) +1;
    nMults++;
  }
  return v[e-1];
}

std::ostream& operator<< (std::ostream &s, const DynamicPtxtPowers &d)
{
  return s << d.getVector();
}

static long
recursivePolyEval(const NTL::ZZX& poly, long k, DynamicPtxtPowers& babyStep,
		  DynamicPtxtPowers& giantStep, long mod,
		  long& recursiveDepth);
static long
PatersonStockmeyer(const NTL::ZZX& poly, long k, long t, long delta,
		   DynamicPtxtPowers& babyStep, DynamicPtxtPowers& giantStep,
		   long mod, long& recursiveDepth);

long simplePolyEval(const NTL::ZZX& poly, DynamicPtxtPowers& babyStep, long mod)
{
  // ensure that we have enough powers
  assertTrue(deg(poly)<=(long)babyStep.size(), "BabyStep has not enough powers (required more than deg(poly))");

  long ret = rem(ConstTerm(poly), mod);
  for (long i=0; i<deg(poly); i++) {
    long coeff = rem(poly[i+1], mod); // f_{i+1}
    long tmp = babyStep[i];           // X^{i+1}
    tmp = MulMod(tmp, coeff, mod);    // f_{i+1} X^{i+1}
    ret = AddMod(ret, tmp, mod);
  }
  return ret;
}

// This procedure assumes that k*(2^e +1) > deg(poly) > k*(2^e -1),
// and that babyStep contains k+ (deg(poly) mod k) powers
static long degPowerOfTwo(const NTL::ZZX& poly, long k, DynamicPtxtPowers& babyStep,
			  DynamicPtxtPowers& giantStep, long mod,
			  long& recursiveDepth)
{
  if (deg(poly)<=babyStep.size()) { // Edge condition, use simple eval
    long ret = simplePolyEval(poly, babyStep, mod);
    recursiveDepth = babyStep.getDepth(deg(poly));
    return ret;
  }
  long subDepth1 =0, subDepth2=0;
  long n = deg(poly)/k;        // We assume n=2^e or n=2^e -1
  n = 1L << NTL::NextPowerOfTwo(n); // round up to n=2^e
  NTL::ZZX r = trunc(poly, (n-1)*k);      // degree <= k(2^e-1)-1
  NTL::ZZX q = RightShift(poly, (n-1)*k); // 0 < degree < 2k
  SetCoeff(r, (n-1)*k);              // monic, degree == k(2^e-1)
  q -= 1;
  if (verbose) std::cerr << ", recursing on "<<r<<" + X^"<<(n-1)*k<<"*"<<q<<std::endl;

  long ret = PatersonStockmeyer(r, k, n/2, 0,
				babyStep, giantStep, mod, subDepth2);
  if (verbose)
    std::cerr << "  PatersonStockmeyer("<<r<<") returns "<<ret
	 << ", depth="<<subDepth2<<std::endl;

  long tmp = simplePolyEval(q, babyStep, mod); // evaluate q
  subDepth1 = babyStep.getDepth(deg(q));
  if (verbose)
    std::cerr << "  simplePolyEval("<<q<<") returns "<<tmp
	 << ", depth="<<subDepth1<<std::endl;

  // multiply by X^{k(n-1)} with minimum depth
  for (long i=1; i<n; i*=2) {
    tmp = MulMod(tmp, giantStep.getPower(i), mod);
    nMults++;
    subDepth1 = max(subDepth1, giantStep.getDepth(i)) +1;
    if (verbose)
      std::cerr << "    after mult by giantStep.getPower("<<i<< ")="
	   << giantStep.getPower(i)<<" of depth="<< giantStep.getDepth(i)
	   << ",  ret="<<tmp<<" and depth is "<<subDepth1<<std::endl;
  }
  totalDepth = max(subDepth1, subDepth2);
  return AddMod(ret, tmp, mod); // return q * X^{k(n-1)} + r
}

// This procedure assumes that poly is monic, deg(poly)=k*(2t-1)+delta
// with t=2^e, and that babyStep contains k+delta powers
static long
PatersonStockmeyer(const NTL::ZZX& poly, long k, long t, long delta,
		   DynamicPtxtPowers& babyStep, DynamicPtxtPowers& giantStep,
		   long mod, long& recursiveDepth)
{
  if (verbose) std::cerr << "PatersonStockmeyer("<<poly<<"), k="<<k<<", t="<<t<<std::endl;

  if (deg(poly)<=babyStep.size()) { // Edge condition, use simple eval
    long ret = simplePolyEval(poly, babyStep, mod);
    recursiveDepth = babyStep.getDepth(deg(poly));
    return ret;
  }
  long subDepth1=0, subDepth2=0;
  long ret, tmp;

  NTL::ZZX r = trunc(poly, k*t);      // degree <= k*2^e-1
  NTL::ZZX q = RightShift(poly, k*t); // degree == k(2^e-1) +delta

  if (verbose) std::cerr << "  r ="<<r<< ", q="<<q;

  const NTL::ZZ& coef = coeff(r,deg(q));
  SetCoeff(r, deg(q), coef-1);  // r' = r - X^{deg(q)}

  if (verbose) std::cerr << ", r'="<<r;

  NTL::ZZX c,s;
  DivRem(c,s,r,q); // r' = c*q + s
  // deg(s)<deg(q), and if c!= 0 then deg(c)<k-delta

  if (verbose) std::cerr << ", c="<<c<< ", s ="<<s<<std::endl;

  assertTrue(deg(s)<deg(q), "Degree of s is not less than degree of q");
  assertTrue(IsZero(c) || deg(c)<k - delta, "Nonzero c has not degree smaller than k - delta");

  SetCoeff(s,deg(q)); // s' = s + X^{deg(q)}, deg(s)==deg(q)

  // reduce the coefficients modulo mod
  for (long i=0; i<=deg(c); i++) rem(c[i],c[i], NTL::to_ZZ(mod));
  c.normalize();
  for (long i=0; i<=deg(s); i++) rem(s[i],s[i], NTL::to_ZZ(mod));
  s.normalize();

  if (verbose) std::cerr << " {t==n+1} recursing on "<<s<<" + (X^"<<t*k<<"+"<<c<<")*"<<q<<std::endl;

  // Evaluate recursively poly = (c+X^{kt})*q + s'
  tmp = simplePolyEval(c, babyStep, mod);
  tmp = AddMod(tmp, giantStep.getPower(t), mod);
  subDepth1 = max(babyStep.getDepth(deg(c)), giantStep.getDepth(t));

  ret = PatersonStockmeyer(q, k, t/2, delta,
			   babyStep, giantStep, mod, subDepth2);

  if (verbose) {
    std::cerr << "  PatersonStockmeyer("<<q<<") returns "<<ret
	 << ", depth="<<subDepth2<<std::endl;
    if (ret != polyEvalMod(q,babyStep[0], mod)) {
      std::stringstream ss;
      ss << "  **1st recursive call failed, q="<<q;
      throw RuntimeError(ss.get());
    }
  }
  ret = MulMod(ret, tmp, mod);
  nMults++;
  subDepth1 = max(subDepth1, subDepth2)+1;

  tmp = PatersonStockmeyer(s, k, t/2, delta,
			   babyStep, giantStep, mod, subDepth2);
  if (verbose) {
    std::cerr << "  PatersonStockmeyer("<<s<<") returns "<<tmp
	 << ", depth="<<subDepth2<<std::endl;
    if (tmp != polyEvalMod(s,babyStep[0], mod)) {
      std::stringstream ss;
      ss << "  **2nd recursive call failed, s="<<s;
      throw RuntimeError(ss.get());
    }
  }
  ret = AddMod(ret,tmp,mod);
  recursiveDepth = max(subDepth1, subDepth2);
  return ret;
}


// This procedure assumes that poly is monic and that babyStep contains
// at least k+delta powers, where delta = deg(poly) mod k
static long
recursivePolyEval(const NTL::ZZX& poly, long k, DynamicPtxtPowers& babyStep,
		  DynamicPtxtPowers& giantStep, long mod,
		  long& recursiveDepth)
{
  if (deg(poly)<=babyStep.size()) { // Edge condition, use simple eval
    long ret = simplePolyEval(poly, babyStep, mod);
    recursiveDepth = babyStep.getDepth(deg(poly));
    return ret;
  }

  if (verbose) std::cerr << "recursivePolyEval("<<poly<<")\n";

  long delta = deg(poly) % k; // deg(poly) mod k
  long n = divc(deg(poly),k); // ceil( deg(poly)/k )
  long t = 1L<<(NTL::NextPowerOfTwo(n)); // t >= n, so t*k >= deg(poly)

  // Special case for deg(poly) = k * 2^e +delta
  if (n==t)
    return degPowerOfTwo(poly, k, babyStep, giantStep, mod, recursiveDepth);

  // When deg(poly) = k*(2^e -1) we use the Paterson-Stockmeyer recursion
  if (n == t-1 && delta==0)
    return PatersonStockmeyer(poly, k, t/2, delta,
			      babyStep, giantStep, mod, recursiveDepth);

  t = t/2;

  // In any other case we have kt < deg(poly) < k(2t-1). We then set
  // u = deg(poly) - k*(t-1) and poly = q*X^u + r with deg(r)<u
  // and recurse on poly = (q-1)*X^u + (X^u+r)

  long u = deg(poly) - k*(t-1);
  NTL::ZZX r = trunc(poly, u);      // degree <= u-1
  NTL::ZZX q = RightShift(poly, u); // degree == k*(t-1)
  q -= 1;
  SetCoeff(r, u);              // degree == u

  long ret, tmp;
  long subDepth1=0, subDepth2=0;
  if (verbose)
    std::cerr << " {deg(poly)="<<deg(poly)<<"<k*(2t-1)="<<k*(2*t-1)
	 << "} recursing on "<<r<<" + X^"<<u<<"*"<<q<<std::endl;
  ret = PatersonStockmeyer(q, k, t/2, 0,
			   babyStep, giantStep, mod, subDepth1);

  if (verbose) {
    std::cerr << "  PatersonStockmeyer("<<q<<") returns "<<ret<<", depth="<<subDepth1<<std::endl;
    if (ret != polyEvalMod(q,babyStep[0], mod)) {
      std::stringstream ss;
      ss << "  @@1st recursive call failed, q="<<q
     	   << ", ret="<<ret<<"!=" << polyEvalMod(q,babyStep[0], mod);
      throw RuntimeError(ss.get());
    }
  }

  tmp = giantStep.getPower(u/k);
  subDepth2 = giantStep.getDepth(u/k);
  if (delta!=0) { // if u is not divisible by k then compute it
    if (verbose)
      std::cerr <<"  multiplying by X^"<<u
	   <<"=giantStep.getPower("<<(u/k)<<")*babyStep.getPower("<<delta<<")="
	   << giantStep.getPower(u/k)<<"*"<<babyStep.getPower(delta)
	   << "="<<tmp<<std::endl;
    tmp = MulMod(tmp, babyStep.getPower(delta), mod);
    nMults++;
    subDepth2++;
  }
  ret = MulMod(ret, tmp, mod);
  nMults ++;
  subDepth1 = max(subDepth1, subDepth2)+1;

  if (verbose) std::cerr << "  after mult by X^{k*"<<u<<"+"<<delta<<"}, depth="<< subDepth1<<std::endl;

  tmp = recursivePolyEval(r, k, babyStep, giantStep, mod, subDepth2);
  if (verbose)
    std::cerr << "  recursivePolyEval("<<r<<") returns "<<tmp<<", depth="<<subDepth2<<std::endl;
  if (tmp != polyEvalMod(r,babyStep[0], mod)) {
    std::stringstream ss;
    ss << "  @@2nd recursive call failed, r="<<r
      << ", ret="<<tmp<<"!=" << polyEvalMod(r,babyStep[0], mod);
    throw RuntimeError(ss.get());
  }
  recursiveDepth = max(subDepth1, subDepth2);
  return AddMod(ret, tmp, mod);
}

// Note: poly is passed by value, not by reference, so the calling routine
// keeps its original polynomial
long evalPolyTopLevel(NTL::ZZX poly, long x, long p, long k=0)
{
  if (verbose)
  std::cerr << "\n* evalPolyTopLevel: p="<<p<<", x="<<x<<", poly="<<poly;

  if (deg(poly)<=2) { // nothing to optimize here
    if (deg(poly)<1) return to_long(coeff(poly, 0));
    DynamicPtxtPowers babyStep(x, p, deg(poly));
    long ret = simplePolyEval(poly, babyStep, p);
    totalDepth = babyStep.getDepth(deg(poly));
    return ret;
  }

  // How many baby steps: set k~sqrt(n/2), rounded up/down to a power of two

  // FIXME: There may be some room for optimization here: it may be possible
  // to choose k as something other than a power of two and still maintain
  // optimal depth, in principle we can try all possible values of k between
  // the two powers of two and choose the one that goves the least number
  // of multiplies, conditioned on minimum depth.

  if (k<=0) {
    long kk = (long) sqrt(deg(poly)/2.0);
    k = 1L << NTL::NextPowerOfTwo(kk);

    // heuristic: if k>>kk then use a smaler power of two
    if ((k==16 && deg(poly)>167) || (k>16 && k>(1.44*kk)))
      k /= 2;
  }
  std::cerr << ", k="<<k;

  long n = divc(deg(poly),k);          // deg(p) = k*n +delta
  if (verbose) std::cerr << ", n="<<n<<std::endl;

  DynamicPtxtPowers babyStep(x, p, k);
  long x2k = babyStep.getPower(k);

  // Special case when deg(p)>k*(2^e -1)
  if (n==(1L << NTL::NextPowerOfTwo(n))) { // n is a power of two
    DynamicPtxtPowers giantStep(x2k, p, n/2, babyStep.getDepth(k));
    if (verbose)
      std::cerr << "babyStep="<<babyStep<<", giantStep="<<giantStep<<std::endl;
    long ret = degPowerOfTwo(poly, k, babyStep, giantStep, p, totalDepth);

    if (verbose) {
      std::cerr << "  degPowerOfTwo("<<poly<<") returns "<<ret<<", depth="<<totalDepth<<std::endl;
      if (ret != polyEvalMod(poly,babyStep[0], p)) {
        std::stringstream ss;
        ss << "  ## recursive call failed, ret="<<ret<<"!="
          << polyEvalMod(poly,babyStep[0], p);
        throw RuntimeError(ss.get());
      }
      // std::cerr << "  babyStep depth=[";
      // for (long i=0; i<babyStep.size(); i++)
      // 	std::cerr << babyStep.getDepth(i+1)<<" ";
      // std::cerr << "]\n";
      // std::cerr << "  giantStep depth=[";
      // for (long i=0; i<giantStep.size(); i++)
      // 	std::cerr<<giantStep.getDepth(i+1)<<" ";
      // std::cerr << "]\n";
    }
    return ret;
  }

  // If n is not a power of two, ensure that poly is monic and that
  // its degree is divisible by k, then call the recursive procedure

  NTL::ZZ topInv; // the inverse mod p of the top coefficient of poly (if any)
  bool divisible = (n*k == deg(poly)); // is the degree divisible by k?
  long nonInvertibe = InvModStatus(topInv, LeadCoeff(poly), NTL::to_ZZ(p));
       // 0 if invertible, 1 if not

  // FIXME: There may be some room for optimization below: instead of
  // adding a term X^{n*k} we can add X^{n'*k} for some n'>n, so long
  // as n' is smaller than the next power of two. We could save a few
  // multiplications since giantStep[n'] may be easier to compute than
  // giantStep[n] when n' has fewer 1's than n in its binary expansion.

  long extra = 0;        // extra!=0 denotes an added term extra*X^{n*k}
  if (!divisible || nonInvertibe) {  // need to add a term
    // set extra = 1 - current-coeff-of-X^{n*k}
    extra = SubMod(1, to_long(coeff(poly,n*k)), p);
    SetCoeff(poly, n*k); // set the top coefficient of X^{n*k} to one
    topInv = NTL::to_ZZ(1);   // inverse of new top coefficient is one
  }

  long t = (extra==0)? divc(n,2) : n;
  DynamicPtxtPowers giantStep(x2k, p, t, babyStep.getDepth(k));

  if (verbose)
    std::cerr << "babyStep="<<babyStep<<", giantStep="<<giantStep<<std::endl;

  long y; // the value to return
  long subDepth1 =0;
  if (!IsOne(topInv)) {
    long top = to_long(poly[n*k]); // record the current top coefficient
    //    std::cerr << ", top-coeff="<<top;

    // Multiply by topInv modulo p to make into a monic polynomial
    poly *= topInv;
    for (long i=0; i<=n*k; i++) rem(poly[i], poly[i], NTL::to_ZZ(p));
    poly.normalize();

    y = recursivePolyEval(poly, k, babyStep, giantStep, p, subDepth1);
    if (verbose) {
      std::cerr << "  recursivePolyEval("<<poly<<") returns "<<y<<", depth="<<subDepth1<<std::endl;
      if (y != polyEvalMod(poly,babyStep[0], p)) {
        std::stringstream ss;
        ss << "## recursive call failed, ret="<<y<<"!="
          << polyEvalMod(poly,babyStep[0], p);
        throw RuntimeError(ss.get());
      }
    }
    y = MulMod(y, top, p); // multiply by the original top coefficient
  }
  else {
    y = recursivePolyEval(poly, k, babyStep, giantStep, p, subDepth1);
    if (verbose) {
      std::cerr << "  recursivePolyEval("<<poly<<") returns "<<y<<", depth="<<subDepth1<<std::endl;
      if (y != polyEvalMod(poly,babyStep[0], p)) {
        std::stringstream ss;
        ss << "## recursive call failed, ret="<<y<<"!="
          << polyEvalMod(poly,babyStep[0], p);
        throw RuntimeError(ss.get());
      }
    }
  }

  if (extra != 0) { // if we added a term, now is the time to subtract back
    if (verbose) std::cerr << ", subtracting "<<extra<<"*X^"<<k*n;
    extra = MulMod(extra, giantStep.getPower(n), p);
    totalDepth = max(subDepth1, giantStep.getDepth(n));
    y = SubMod(y, extra, p);
  }
  else totalDepth = subDepth1;
  if (verbose) std::cerr << std::endl;
  return y;
}
#endif


// ============================================================
// Asymmetric BSGS Polynomial Evaluation (Wang & Wang, ASIA CCS 2025)
// ============================================================

static bool ispowof2(long i) { return i > 0 && (i & (i - 1)) == 0; }
static long prepowof2(long i) { return 1l << long(std::floor(std::log2(i))); }

struct CtxtMap
{
  std::map<long, long> idx;     // std::map 要求 value 存在隐式构造器
  std::vector<Ctxt> vec; // 使用 vector 模拟 map

  void insert(long i, const Ctxt& ctxt)
  {
    idx[i] = vec.size();
    vec.push_back(ctxt);
  }

  bool find(long i) { return idx.find(i) != idx.end(); }

  Ctxt& operator[](long i) { return vec[idx[i]]; }
};

struct Box
{
  long factor;
  long boxsize; // idx = [factor, 2*factor, ..., boxsize*factor]

  Box(long factor, long boxsize) : factor(factor), boxsize{boxsize} {}
};

struct Boxes
{
  long deg;
  long logdeg;
  long num;
  std::vector<Box> boxes;
  CtxtMap map;

  Boxes(long deg, long num) : deg(deg), num(num)
  {
    boxes.clear();
    logdeg = ceil(log2(NTL::xdouble(deg)));
    long div = logdeg / num; // 基础长度，均匀分配各个 box 大小
    long rem = logdeg % num; // 前 rem <= num-1 个盒子加长 1 比特
    long len = 0;

    long j = 0;
    for (; j < rem; j++) {
      boxes.push_back(Box(1l << len, (1l << (div + 1)) - 1)); // 前 rem 个盒子
      len += div + 1;
    }
    for (; j < num - 1; j++) {
      boxes.push_back(Box(1l << len, (1l << div) - 1)); // 中间 num-rem-1 个盒子
      len += div;
    }
    boxes.push_back(Box(1l << len, deg >> len)); // 最后一个盒子
  }

  void insert(long i, const Ctxt& ctxt) { map.insert(i, ctxt); }

  bool find(long i) { return map.find(i); }

  Ctxt& operator[](long i) { return map[i]; }

  void decomp(long i, std::vector<long>& idx)
  {
    idx.clear();
    for (long k = 0; k < num; k++) {
      long j =
          (i / boxes[k].factor) % (boxes[k].boxsize + 1); // 分解为 Boxes 的索引
      if (j != 0)
        idx.push_back(j * boxes[k].factor);
    }
  }
};

IndexSet findcommonPrimeSet(
    Ctxt& ctxt) // 计算张量积之前，控制噪声
{
  const double safety = 1 * log(2.0);
  const double slack = 4 * log(2.0);

  double cap = ctxt.logOfPrimeSet() -
               NTL::log(std::max(ctxt.getNoiseBound(), NTL::to_xdouble(1.0)));
  double adn = log(ctxt.modSwitchAddedNoiseBound());

  double lo, hi;
  hi = cap + adn - safety;
  lo = hi - slack;

  return ctxt.context.getModSizeTable().getSet4Size(lo, hi, ctxt.primeSet, 0);
}

void mulAtSameLevel(Ctxt &res, Ctxt &ct1, Ctxt &ct2, IndexSet &Primes)
{
    ct1.bringToSet(Primes);
    ct2.bringToSet(Primes);

    Ctxt tmp = ct1;
    res.tensorProduct(tmp, ct2);
}

void BS(Ctxt& res, const std::vector<Ctxt>& PowsBS, NTL::ZZX f)
{
  long deg = NTL::deg(f);
  if (deg == -1) // 零多项式
  {
    return;
  }

  res.addConstant(f[0]);
  for (long i = 1; i <= deg; i++) {
    if (f[i] != 0) {
      Ctxt tmp =
          PowsBS[i -
                 1]; // 假设 powsBS 已经计算所必须的幂次，x^i 存储在 [i-1] 位置
      tmp.multByConstant(f[i]);
      res.addCtxt(tmp);
    }
  }
}

void GS(Ctxt& res,
        const std::vector<Ctxt>& PowsBS,
        const std::vector<Ctxt>& PowsGS,
        long k,
        const Ctxt& empty,
        NTL::ZZX f)
{
  long deg = NTL::deg(f);
  if (deg == -1) // 零多项式
  {
    return;
  } else if (deg == 0) // 常数多项式
  {
    res.addConstant(f[0]);
    return;
  }

  long n = deg / k; // floor(deg/k), n*k <= deg
  Ctxt subf(empty);
  BS(res, PowsBS, NTL::trunc(f, k));
  f = NTL::RightShift(f, k);

  for (long j = 1; j < n; j++) {
    Ctxt tmp(empty);
    BS(tmp, PowsBS, NTL::trunc(f, k)); // [j*k, ..., (j+1)*k-1], j >= 1
    f = NTL::RightShift(f, k);
    subf.tensorProduct(tmp, PowsGS[j - 1]); // x^{j*k}
    res.addCtxt(subf);
  }

  if (deg >= n * k) {
    Ctxt tmp(empty);
    BS(tmp, PowsBS, f);                     // [n*k, ..., deg]
    subf.tensorProduct(tmp, PowsGS[n - 1]); // x^{n*k}
    res.addCtxt(subf);
  }
}

long MakeBSGS(long& degreeHigh,
              long& degreeLow,
              long& k,
              CtxtMap& powHigh,
              std::vector<Ctxt>& PowsBS,
              std::vector<Ctxt>& PowsGS,
              Ctxt& empty,
              const Ctxt& x,
              long deg,
              long num)
{
  long logdeg = ceil(log2(NTL::xdouble(
      deg))); // 注意不要用 log2(deg+1)，因为 2^i 的乘法深度依旧是 i 而非 i+1
  while (logdeg < 2 * num - 1) // 本算法的最小适用范围
    num -= 1;

  if (num == 1)
    return num;

  long factorHigh = 1l << long(logdeg - (num - 1)); // 最高 num-1 比特的间距
  degreeHigh = deg >> long(logdeg - (num - 1)); // 最高 num-1 比特的范围
  degreeLow = factorHigh - 1; // 较低的比特做 BSGS 不做模切换的版本

  Boxes BOXES(degreeLow, num); // 建立 num 个集合，直和为 [0, degreeLow]
  // CtxtMap powHigh;             // factorHigh 的倍数, [1, degreeHigh]

  for (long k = 0; k < num; k++) {
    long fact = BOXES.boxes[k].factor;
    long size = BOXES.boxes[k].boxsize;

    for (long j = 1; j <= size; j++) {
      long i = j * fact; // 仅计算这些索引的 x^i
      if (i == 1) {
        BOXES.insert(1, x); // 存储 x^1
        continue;
      }

      if (ispowof2(i)) {
        BOXES.insert(i, BOXES[i >> 1]); // 幂次是二的幂次
        BOXES[i].square();
      } else {
        long powof2 = prepowof2(i);
        BOXES.insert(i, BOXES[powof2]);
        BOXES[i].multiplyBy(BOXES[i - powof2]);
      }
    }
  }

  powHigh.insert(1, BOXES[factorHigh >> 1]); // x^factorHigh
  powHigh[1].square();

  for (long i = 2; i <= degreeHigh;
       i++) // [factorHigh, 2*factorHigh, ..., degreeHigh*factorHigh]
  {
    if (ispowof2(i)) {
      powHigh.insert(i, powHigh[i >> 1]);
      powHigh[i].square();
    } else {
      long powof2 = prepowof2(i);
      powHigh.insert(i, powHigh[powof2]);
      powHigh[i].multiplyBy(powHigh[i - powof2]);
    }
  }

  IndexSet Primes =
      findcommonPrimeSet(BOXES.map.vec.back()); // BOXES 中的最深处密文
  // Ctxt empty(ZeroCtxtLike, x);                         //
  // 空密文的模板
  empty.clear();
  empty.primeSet = Primes;

  for (Ctxt& ctxt : BOXES.map.vec) // 所有密文调整到同一层
  {
    ctxt.bringToSet(Primes);
  }

  k = BOXES.boxes[BOXES.num / 2].factor; // 前 num/2 个盒子的范围
  // std::vector<Ctxt> PowsBS(k - 1, empty);         // [1, 2, ..., k-1]
  // std::vector<Ctxt> PowsGS(degreeLow / k, empty); // [k, 2k, ...,
  // (degreeLow/k)*k]
  PowsBS.resize(k - 1, empty);
  PowsGS.resize(degreeLow / k, empty);

  for (long i = 1; i < k; i++) {
    std::vector<long> idx; // 把 i 分解
    BOXES.decomp(i, idx);

    long l = idx.size();
    Ctxt prod = BOXES[idx[0]]; // 集合 idx 非空
    for (long j = 1; j < l; j++) {
      Ctxt tmp(empty);
      tmp.tensorProduct(prod, BOXES[idx[j]]);
      prod = tmp;
    }

    PowsBS[i - 1] = prod; // x^i
  }

  for (long i = k; i <= degreeLow; i += k) {
    std::vector<long> idx; // 把 i 分解
    BOXES.decomp(i, idx);

    long l = idx.size();
    Ctxt prod = BOXES[idx[0]]; // 集合 idx 非空
    for (long j = 1; j < l; j++) {
      Ctxt tmp(empty);
      tmp.tensorProduct(prod, BOXES[idx[j]]);
      prod = tmp;
    }

    PowsGS[i / k - 1] = prod; // x^i
  }

  return num;
}

void AsymBSGB(Ctxt& res,
              Ctxt x,
              NTL::ZZX f,
              long& degreeHigh,
              long& degreeLow,
              long& k,
              CtxtMap& powHigh,
              std::vector<Ctxt>& PowsBS,
              std::vector<Ctxt>& PowsGS,
              Ctxt& empty,
              long num)
{
    x.cleanUp();

    long deg = NTL::deg(f);
    if (num <= 1) // 经典算法
    {
        std::vector<Ctxt> pows;
        pows.push_back(x);
        for (long i = 2; i <= deg; i++)
        {
            if (ispowof2(i))
            {
                pows.push_back(pows[(i >> 1) - 1]);
                pows[i - 1].square();
            }
            else
            {
                long powof2 = prepowof2(i);
                pows.push_back(pows[powof2 - 1]);
                pows[i - 1].multiplyBy(pows[i - powof2 - 1]);
            }
        }

        res.clear();
        BS(res, pows, f);
        return;
    }


    res = empty;
    GS(res, PowsBS, PowsGS, k, empty, NTL::trunc(f, degreeLow + 1));
    f = NTL::RightShift(f, degreeLow + 1);


    // IndexSet level = findcommonPrimeSet(res);
    // res.modDownToSet(level);

    // cout << res.primeSet << endl;
    // printf("bit Capacity of res2 = %ld\n", res.bitCapacity());
    // where;

    for (long i = 1; i < degreeHigh; i++)
    {
        Ctxt subf(empty);
        GS(subf, PowsBS, PowsGS, k, empty, NTL::trunc(f, degreeLow + 1));
        f = NTL::RightShift(f, degreeLow + 1);

        subf.multLowLvl(powHigh[i]); // 自动模切换，不做重线性化
        //mulAtSameLevel(subf, subf, powHigh[i], level);

        res.modDownToSet(subf.primeSet);
        res.addCtxt(subf);
    }

    Ctxt subf(empty);
    GS(subf, PowsBS, PowsGS, k, empty, f);

    subf.multLowLvl(powHigh[degreeHigh]); // 自动模切换，不做重线性化
    //mulAtSameLevel(subf, subf, powHigh[degreeHigh], level);

    res.modDownToSet(subf.primeSet);
    res.addCtxt(subf);

    res.bringToSet(x.primeSet);
    
    res.reLinearize();
    // res.bringToSet(level);
    // cout << res.primeSet << endl;
    // printf("bit Capacity of res = %ld\n", res.bitCapacity());
    // where;
}

void PolyEvalAsymBSGB(Ctxt& res,
                      const Ctxt& x,
                      NTL::ZZX f,
                      long num)
{
    long deg = NTL::deg(f);
    if (deg <= 1)
    {
        if(deg == 1){
            res = x;
            res.multByConstant(f[1]);
        }
        else
            res.clear(); // 在 deg=1 时直接复制，不必执行 clear

        res.addConstant(f[0]);
        return;
    }

    CtxtMap powHigh;
    std::vector<Ctxt> PowsBS, PowsGS;
    Ctxt empty(ZeroCtxtLike, x);
    long degreeHigh, degreeLow, k;

    num = MakeBSGS(degreeHigh, degreeLow, k, powHigh, PowsBS, PowsGS, empty, x, deg, num);
    AsymBSGB(res, x, f, degreeHigh, degreeLow, k, powHigh, PowsBS, PowsGS, empty, num);

  //cout << x.primeSet << endl;
  //cout << res.primeSet << endl;
}

void ManyPolyEvalAsymBSGB(std::vector<Ctxt*>& res,
                          const Ctxt& x,
                          NTL::vec_ZZX fs,
                          long num)
{
  long len = fs.length();
    long max_degree = 0;

    for (auto &f : fs)
    {
        long degree = NTL::deg(f);
        max_degree = (degree > max_degree) ? degree : max_degree;
    }

    CtxtMap powHigh;
    std::vector<Ctxt> PowsBS, PowsGS;
    Ctxt empty(ZeroCtxtLike, x);
    long degreeHigh, degreeLow, k;

    num = MakeBSGS(degreeHigh, degreeLow, k, powHigh, PowsBS, PowsGS, empty, x, max_degree, num);

    for (long i = 0; i < len; i++)
    {
        long deg = NTL::deg(fs[i]);
        if (deg == -1)
        {
            res.clear();
            return;
        }

        // NTL::rem(fs[i], fs[i], NTL::ZZX(x.getPtxtSpace()));

        if (deg <= 1)
        {

            if (deg == 1)
            {
                *res[i] = x;
                res[i]->multByConstant(fs[i][1]);
            }
            else
                res.clear(); // 在 deg=1 时直接复制，不必执行 clear

            res[i]->addConstant(fs[i][0]);
            return;
        }

        AsymBSGB(*res[i], x, fs[i], degreeHigh, degreeLow, k, powHigh, PowsBS, PowsGS, empty, num);
    }
}

// ============================================================

} // namespace helib
