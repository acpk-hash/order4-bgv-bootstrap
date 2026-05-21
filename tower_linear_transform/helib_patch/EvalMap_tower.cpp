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
#include <helib/EvalMap.h>
#include <helib/apiAttributes.h>

// needed to get NTL's TraceMap functions...needed for ThinEvalMap
#include <NTL/lzz_pXFactoring.h>
#include <NTL/GF2XFactoring.h>

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-local-typedefs"
#pragma GCC diagnostic ignored "-Wunused-variable"
#pragma GCC diagnostic ignored "-Wunused-parameter"

namespace helib {

// Forward declarations
static void towerBenchmark96(const EncryptedArray& ea,
                             const NTL::Vec<long>& reps,
                             long dim,
                             long cofactor);
static BlockMatMul1D* buildStep1Matrix(const EncryptedArray& ea,
                                       std::shared_ptr<CubeSignature> sig,
                                       const NTL::Vec<long>& reps,
                                       long dim,
                                       long cofactor,
                                       bool invert,
                                       bool normal_basis);
static MatMul1D* buildStep2Matrix(const EncryptedArray& ea,
                                  std::shared_ptr<CubeSignature> sig,
                                  const NTL::Vec<long>& reps,
                                  long dim,
                                  long cofactor,
                                  bool invert);
static void init_representatives(NTL::Vec<long>& representatives,
                                 long dim,
                                 const NTL::Vec<long>& mvec,
                                 const PAlgebra& zMStar);

// Constructor: initializing tables for the evaluation-map transformations

EvalMap::EvalMap(const EncryptedArray& _ea,
                 bool minimal,
                 const NTL::Vec<long>& mvec,
                 bool _invert,
                 bool build_cache,
                 bool normal_basis) :
    ea(_ea), invert(_invert)
{
  const PAlgebra& zMStar = ea.getPAlgebra();

  long p = zMStar.getP();
  long d = zMStar.getOrdP();

  // FIXME: we should check that ea was initialized with
  // G == factors[0], but this is a slight pain to check
  // currently

  // NOTE: this code is derived from a more general setting, and
  // could certainly be greatly simplified

  nfactors = mvec.length();
  assertTrue(nfactors > 0, "Invalid argument: mvec must not be empty");

  for (long i = 0; i < nfactors; i++) {
    for (long j = i + 1; j < nfactors; j++) {
      assertEq(NTL::GCD(mvec[i], mvec[j]),
               1l,
               "Invalid argument: mvec elements must be pairwise co-prime");
    }
  }

  long m = computeProd(mvec);
  assertEq(m,
           (long)zMStar.getM(),
           "Invalid argument: Product of mvec elements does not match "
           "ea.zMStar.getM()");

  NTL::Vec<long> phivec(NTL::INIT_SIZE, nfactors);
  for (long i = 0; i < nfactors; i++)
    phivec[i] = phi_N(mvec[i]);
  long phim = computeProd(phivec);

  NTL::Vec<long> dprodvec(NTL::INIT_SIZE, nfactors + 1);
  dprodvec[nfactors] = 1;

  for (long i = nfactors - 1; i >= 0; i--)
    dprodvec[i] =
        dprodvec[i + 1] *
        multOrd(NTL::PowerMod(p % mvec[i], dprodvec[i + 1], mvec[i]), mvec[i]);

  NTL::Vec<long> dvec(NTL::INIT_SIZE, nfactors);
  for (long i = 0; i < nfactors; i++)
    dvec[i] = dprodvec[i] / dprodvec[i + 1];

  long nslots = phim / d;
  assertEq(d, dprodvec[0], "dprodvec must start with d");
  assertEq(nslots,
           (long)zMStar.getNSlots(),
           "Slot count mismatch between ea and phi(m)/d");

  long inertPrefix = 0;
  for (long i = 0; i < nfactors && dvec[i] == 1; i++) {
    inertPrefix++;
  }

  if (inertPrefix != nfactors - 1)
    throw LogicError("EvalMap: case not handled: bad inertPrefix");

  NTL::Vec<NTL::Vec<long>> local_reps(NTL::INIT_SIZE, nfactors);
  for (long i = 0; i < nfactors; i++)
    init_representatives(local_reps[i], i, mvec, zMStar);

  NTL::Vec<long> crtvec(NTL::INIT_SIZE, nfactors);
  for (long i = 0; i < nfactors; i++)
    crtvec[i] = (m / mvec[i]) * NTL::InvMod((m / mvec[i]) % mvec[i], mvec[i]);

  NTL::Vec<long> redphivec(NTL::INIT_SIZE, nfactors);
  for (long i = 0; i < nfactors; i++)
    redphivec[i] = phivec[i] / dvec[i];

  CubeSignature redphisig(redphivec);

  NTL::Vec<std::shared_ptr<CubeSignature>> sig_sequence;
  sig_sequence.SetLength(nfactors + 1);
  sig_sequence[nfactors] = std::make_shared<CubeSignature>(phivec);

  NTL::Vec<long> reduced_phivec = phivec;

  for (long dim = nfactors - 1; dim >= 0; dim--) {
    reduced_phivec[dim] /= dvec[dim];
    sig_sequence[dim] = std::make_shared<CubeSignature>(reduced_phivec);
  }

  long dim = nfactors - 1;
  std::unique_ptr<BlockMatMul1D> mat1_data;
  mat1_data.reset(buildStep1Matrix(ea,
                                   sig_sequence[dim],
                                   local_reps[dim],
                                   dim,
                                   m / mvec[dim],
                                   invert,
                                   normal_basis));
  mat1.reset(new BlockMatMul1DExec(*mat1_data, minimal));

  matvec.SetLength(nfactors - 1);
  for (dim = nfactors - 2; dim >= 0; --dim) {
    std::unique_ptr<MatMul1D> mat_data;

    mat_data.reset(buildStep2Matrix(ea,
                                    sig_sequence[dim],
                                    local_reps[dim],
                                    dim,
                                    m / mvec[dim],
                                    invert));
    matvec[dim].reset(new MatMul1DExec(*mat_data, minimal));
  }

  if (build_cache)
    upgrade();
}

void EvalMap::upgrade()
{
  mat1->upgrade();
  for (long i = 0; i < matvec.length(); i++)
    matvec[i]->upgrade();
}

// Applying the evaluation (or its inverse) map to a ciphertext
void EvalMap::apply(Ctxt& ctxt) const
{
  if (!invert) { // forward direction
    mat1->mul(ctxt);

    for (long i = matvec.length() - 1; i >= 0; i--)
      matvec[i]->mul(ctxt);
  } else { // inverse transformation
    for (long i = 0; i < matvec.length(); i++)
      matvec[i]->mul(ctxt);

    mat1->mul(ctxt);
  }
}

static void init_representatives(NTL::Vec<long>& representatives,
                                 long dim,
                                 const NTL::Vec<long>& mvec,
                                 const PAlgebra& zMStar)
{
  assertInRange(dim,
                0l,
                mvec.length(),
                "Invalid argument: dim must be between 0 and mvec.length()");

  // special case
  if (dim >= LONG(zMStar.numOfGens())) {
    representatives.SetLength(1);
    representatives[0] = 1;
    return;
  }

  long m = mvec[dim];
  long D = zMStar.OrderOf(dim);
  long g = NTL::InvMod(zMStar.ZmStarGen(dim) % m, m);

  representatives.SetLength(D);
  for (long i = 0; i < D; i++)
    representatives[i] = NTL::PowerMod(g, i, m);
}

// The callback interface for the matrix-multiplication routines.

//! \cond FALSE (make doxygen ignore these classes)
template <typename type>
class Step2Matrix : public MatMul1D_derived<type>
{
  PA_INJECT(type)

  const EncryptedArray& base_ea;
  std::shared_ptr<CubeSignature> sig;
  long dim;
  NTL::Mat<RX> A;

public:
  // constructor
  Step2Matrix(const EncryptedArray& _ea,
              std::shared_ptr<CubeSignature> _sig,
              const NTL::Vec<long>& reps,
              long _dim,
              long cofactor,
              bool invert = false) :
      base_ea(_ea), sig(_sig), dim(_dim)
  {
    long sz = sig->getDim(dim);
    assertEq(sz,
             reps.length(),
             "Invalid argument: sig->getDim(dim) must equal reps.length()");

    const EncryptedArrayDerived<type>& ea = _ea.getDerived(type());
    RBak bak;
    bak.save();
    _ea.getAlMod().restoreContext();
    const RX& G = ea.getG();

    NTL::Vec<RX> points(NTL::INIT_SIZE, sz);
    for (long j = 0; j < sz; j++)
      points[j] = RX(reps[j] * cofactor, 1) % G;

    A.SetDims(sz, sz);
    for (long j = 0; j < sz; j++)
      A[0][j] = 1;

    for (long i = 1; i < sz; i++)
      for (long j = 0; j < sz; j++)
        A[i][j] = (A[i - 1][j] * points[j]) % G;

    if (invert) {
      REBak ebak;
      ebak.save();
      ea.restoreContextForG();

      mat_RE A1, A2;
      conv(A1, A);

      long p = _ea.getAlMod().getZMStar().getP();
      long r = _ea.getAlMod().getR();

      ppInvert(A2, A1, p, r);
      conv(A, A2);
    }
  }

  bool get(RX& out, long i, long j, UNUSED long k) const override
  {
    out = A[i][j];
    return false;
  }

  const EncryptedArray& getEA() const override { return base_ea; }
  bool multipleTransforms() const override { return false; }
  long getDim() const override { return dim; }
};

static MatMul1D* buildStep2Matrix(const EncryptedArray& ea,
                                  std::shared_ptr<CubeSignature> sig,
                                  const NTL::Vec<long>& reps,
                                  long dim,
                                  long cofactor,
                                  bool invert)
{
  switch (ea.getTag()) {
  case PA_GF2_tag:
    return new Step2Matrix<PA_GF2>(ea, sig, reps, dim, cofactor, invert);

  case PA_zz_p_tag:
    return new Step2Matrix<PA_zz_p>(ea, sig, reps, dim, cofactor, invert);

  default:
    return 0;
  }
}

template <typename type>
class Step1Matrix : public BlockMatMul1D_derived<type>
{
  PA_INJECT(type)

  const EncryptedArray& base_ea;
  std::shared_ptr<CubeSignature> sig;
  long dim;
  NTL::Mat<mat_R> A;

public:
  // constructor
  Step1Matrix(const EncryptedArray& _ea,
              std::shared_ptr<CubeSignature> _sig,
              const NTL::Vec<long>& reps,
              long _dim,
              long cofactor,
              bool invert,
              bool normal_basis) :
      base_ea(_ea), sig(_sig), dim(_dim)
  {
    const EncryptedArrayDerived<type>& ea = _ea.getDerived(type());
    RBak bak;
    bak.save();
    _ea.getAlMod().restoreContext();
    const RX& G = ea.getG();
    long d = deg(G);

    long sz = sig->getDim(dim);
    assertEq(sz,
             reps.length(),
             "Invalid argument: sig->getDim(dim) must equal reps.length()");
    assertEq(dim,
             sig->getNumDims() - 1,
             "Invalid argument: dim must be one less than sig->getNumDims()");
    assertEq(sig->getSize(), ea.size(), "sig and ea do not have matching size");

    // so sz == phi(m_last)/d, where d = deg(G) = order of p mod m

    NTL::Vec<RX> points(NTL::INIT_SIZE, sz);
    for (long j = 0; j < sz; j++)
      points[j] = RX(reps[j] * cofactor, 1) % G;

    NTL::Mat<RX> AA(NTL::INIT_SIZE, sz * d, sz);
    for (long j = 0; j < sz; j++)
      AA[0][j] = 1;

    for (long i = 1; i < sz * d; i++)
      for (long j = 0; j < sz; j++)
        AA[i][j] = (AA[i - 1][j] * points[j]) % G;

    A.SetDims(sz, sz);
    for (long i = 0; i < sz; i++)
      for (long j = 0; j < sz; j++) {
        A[i][j].SetDims(d, d);
        for (long k = 0; k < d; k++)
          VectorCopy(A[i][j][k], AA[i * d + k][j], d);
      }

    if (invert) {
      mat_R A1, A2;
      A1.SetDims(sz * d, sz * d);
      for (long i = 0; i < sz * d; i++)
        for (long j = 0; j < sz * d; j++)
          A1[i][j] = A[i / d][j / d][i % d][j % d];

      long p = _ea.getAlMod().getZMStar().getP();
      long r = _ea.getAlMod().getR();

      ppInvert(A2, A1, p, r);

      for (long i = 0; i < sz * d; i++)
        for (long j = 0; j < sz * d; j++)
          A[i / d][j / d][i % d][j % d] = A2[i][j];

      if (normal_basis) {
        const NTL::Mat<R>& CB = ea.getNormalBasisMatrix();

        // multiply each entry of A on the right by CB
        for (long i = 0; i < sz; i++)
          for (long j = 0; j < sz; j++)
            A[i][j] = A[i][j] * CB;
      } // if (normal_basis)
    }   // if (invert)
  }     // constructor

  bool get(mat_R& out, long i, long j, UNUSED long k) const override
  {
    out = A[i][j];
    return false;
  }

  const EncryptedArray& getEA() const override { return base_ea; }
  bool multipleTransforms() const override { return false; }
  long getDim() const override { return dim; }
};

static BlockMatMul1D* buildStep1Matrix(const EncryptedArray& ea,
                                       std::shared_ptr<CubeSignature> sig,
                                       const NTL::Vec<long>& reps,
                                       long dim,
                                       long cofactor,
                                       bool invert,
                                       bool normal_basis)
{
  switch (ea.getTag()) {
  case PA_GF2_tag:
    return new Step1Matrix<PA_GF2>(ea,
                                   sig,
                                   reps,
                                   dim,
                                   cofactor,
                                   invert,
                                   normal_basis);

  case PA_zz_p_tag:
    return new Step1Matrix<PA_zz_p>(ea,
                                    sig,
                                    reps,
                                    dim,
                                    cofactor,
                                    invert,
                                    normal_basis);

  default:
    return 0;
  }
}
//! \endcond

//=============== ThinEvalMap stuff

// needed to make generic programming work

void RelaxedInv(NTL::Mat<NTL::zz_p>& x, const NTL::Mat<NTL::zz_p>& a)
{
  relaxed_inv(x, a);
}

void RelaxedInv(NTL::Mat<NTL::GF2>& x, const NTL::Mat<NTL::GF2>& a)
{
  inv(x, a);
}

void TraceMap(NTL::GF2X& w,
              const NTL::GF2X& a,
              long d,
              const NTL::GF2XModulus& F,
              const NTL::GF2X& b)

{
  if (d < 0)
    throw InvalidArgument("TraceMap: d is negative");

  NTL::GF2X y, z, t;

  z = b;
  y = a;
  clear(w);

  while (d) {
    if (d == 1) {
      if (IsZero(w))
        w = y;
      else {
        CompMod(w, w, z, F);
        add(w, w, y);
      }
    } else if ((d & 1) == 0) {
      Comp2Mod(z, t, z, y, z, F);
      add(y, t, y);
    } else if (IsZero(w)) {
      w = y;
      Comp2Mod(z, t, z, y, z, F);
      add(y, t, y);
    } else {
      Comp3Mod(z, t, w, z, y, w, z, F);
      add(w, w, y);
      add(y, t, y);
    }

    d = d >> 1;
  }
}

// Forward declarations
static MatMul1D* buildThinStep1Matrix(const EncryptedArray& ea,
                                      std::shared_ptr<CubeSignature> sig,
                                      const NTL::Vec<long>& reps,
                                      long dim,
                                      long cofactor);
static MatMul1D* buildThinStep2Matrix(const EncryptedArray& ea,
                                      std::shared_ptr<CubeSignature> sig,
                                      const NTL::Vec<long>& reps,
                                      long dim,
                                      long cofactor,
                                      bool invert,
                                      bool inflate = false);
static void init_representatives(NTL::Vec<long>& representatives,
                                 long dim,
                                 const NTL::Vec<long>& mvec,
                                 const PAlgebra& zMStar);

// Tower 2-Stage class (defined early for use in constructor)
template <typename type>
class ThinStep2Tower96Stage : public MatMul1D_derived<type>
{
  PA_INJECT(type)
  const EncryptedArray& base_ea;
  long dim;
  NTL::Mat<RX> A;
public:
  ThinStep2Tower96Stage(const EncryptedArray& _ea, long _dim,
                        const NTL::Mat<RX>& _A)
      : base_ea(_ea), dim(_dim), A(_A) {}

  bool get(RX& out, long i, long j, UNUSED long k) const override {
    out = A[i][j];
    return IsZero(out);
  }
  const EncryptedArray& getEA() const override { return base_ea; }
  bool multipleTransforms() const override { return false; }
  long getDim() const override { return dim; }
};

template <typename type>
static void buildTowerStages96(const EncryptedArray& _ea,
                               const NTL::Vec<long>& reps,
                               long dim,
                               long cofactor,
                               NTL::Mat<typename type::RX>& M1_out,
                               NTL::Mat<typename type::RX>& M2_out,
                               long& m2_ndiags);

template <typename type>
static void buildButterflyStages96(const EncryptedArray& _ea,
                                    const NTL::Vec<long>& reps,
                                    long dim,
                                    long cofactor,
                                    bool invert,
                                    std::vector<NTL::Mat<typename type::RX>>& stages_out);

template <typename type>
static void buildRaderVandermonde96(const EncryptedArray& _ea,
                                    const NTL::Vec<long>& reps,
                                    long dim, long cofactor, bool invert,
                                    NTL::Mat<typename type::RX>& V_out);

template <typename type>
static void buildRaderFactored96(const EncryptedArray& _ea,
                                  const NTL::Vec<long>& reps,
                                  long dim, long cofactor, bool invert,
                                  std::vector<NTL::Mat<typename type::RX>>& stages_out);

// Constructor: initializing tables for the evaluation-map transformations

ThinEvalMap::ThinEvalMap(const EncryptedArray& _ea,
                         bool minimal,
                         const NTL::Vec<long>& mvec,
                         bool _invert,
                         bool build_cache) :
    ea(_ea), invert(_invert)
{
  const PAlgebra& zMStar = ea.getPAlgebra();

  long p = zMStar.getP();
  long d = zMStar.getOrdP();
  long sz = zMStar.numOfGens();

  // FIXME: we should check that ea was initialized with
  // G == factors[0], but this is a slight pain to check
  // currently

  // NOTE: this code is derived from a more general setting, and
  // could certainly be greatly simplified

  nfactors = mvec.length();
  assertTrue(nfactors > 0, "Invalid argument: mvec must have positive length");

  for (long i = 0; i < nfactors; i++) {
    for (long j = i + 1; j < nfactors; j++) {
      assertEq(NTL::GCD(mvec[i], mvec[j]),
               1l,
               "Invalid argument: mvec must have pairwise-disjoint entries");
    }
  }

  long m = computeProd(mvec);
  assertEq(m,
           (long)zMStar.getM(),
           "Invalid argument: mvec's product does not match ea's m");

  NTL::Vec<long> phivec(NTL::INIT_SIZE, nfactors);
  for (long i = 0; i < nfactors; i++)
    phivec[i] = phi_N(mvec[i]);
  long phim = computeProd(phivec);

  NTL::Vec<long> dprodvec(NTL::INIT_SIZE, nfactors + 1);
  dprodvec[nfactors] = 1;

  for (long i = nfactors - 1; i >= 0; i--)
    dprodvec[i] =
        dprodvec[i + 1] *
        multOrd(NTL::PowerMod(p % mvec[i], dprodvec[i + 1], mvec[i]), mvec[i]);

  NTL::Vec<long> dvec(NTL::INIT_SIZE, nfactors);
  for (long i = 0; i < nfactors; i++)
    dvec[i] = dprodvec[i] / dprodvec[i + 1];

  long nslots = phim / d;
  assertEq(d, dprodvec[0], "d must match the first entry of dprodvec");
  assertEq(nslots,
           (long)zMStar.getNSlots(),
           "Invalid argument: mismatch of number of slots");

  long inertPrefix = 0;
  for (long i = 0; i < nfactors && dvec[i] == 1; i++) {
    inertPrefix++;
  }

  if (inertPrefix != nfactors - 1)
    throw LogicError("ThinEvalMap: case not handled: bad inertPrefix");

  NTL::Vec<NTL::Vec<long>> local_reps(NTL::INIT_SIZE, nfactors);
  for (long i = 0; i < nfactors; i++)
    init_representatives(local_reps[i], i, mvec, zMStar);

  NTL::Vec<long> crtvec(NTL::INIT_SIZE, nfactors);
  for (long i = 0; i < nfactors; i++)
    crtvec[i] = (m / mvec[i]) * NTL::InvMod((m / mvec[i]) % mvec[i], mvec[i]);

  NTL::Vec<long> redphivec(NTL::INIT_SIZE, nfactors);
  for (long i = 0; i < nfactors; i++)
    redphivec[i] = phivec[i] / dvec[i];

  CubeSignature redphisig(redphivec);

  NTL::Vec<std::shared_ptr<CubeSignature>> sig_sequence;
  sig_sequence.SetLength(nfactors + 1);
  sig_sequence[nfactors] = std::make_shared<CubeSignature>(phivec);

  NTL::Vec<long> reduced_phivec = phivec;

  for (long dim = nfactors - 1; dim >= 0; dim--) {
    reduced_phivec[dim] /= dvec[dim];
    sig_sequence[dim] = std::make_shared<CubeSignature>(reduced_phivec);
  }

  matvec.SetLength(nfactors);

  std::cerr << "[ThinEvalMap] nfactors=" << nfactors << ", sz=" << sz
            << ", invert=" << invert << "\n";
  for (long i = 0; i < nfactors; i++)
    std::cerr << "  dim " << i << ": local_reps.length()=" << local_reps[i].length()
              << ", mvec[i]=" << mvec[i] << "\n";

  if (invert) {
    long dim = nfactors - 1;
    std::unique_ptr<MatMul1D> mat1_data;
    mat1_data.reset(buildThinStep1Matrix(ea,
                                         sig_sequence[dim],
                                         local_reps[dim],
                                         dim,
                                         m / mvec[dim]));
    matvec[dim].reset(new MatMul1DExec(*mat1_data, minimal));
  } else if (sz == nfactors) {
    long dim = nfactors - 1;
    std::cerr << "[ThinEvalMap] sz==nfactors branch: dim=" << dim
              << ", local_reps[dim].length()=" << local_reps[dim].length()
              << ", mvec[dim]=" << mvec[dim]
              << ", m/mvec[dim]=" << m/mvec[dim]
              << ", invert=" << invert << "\n";
    std::unique_ptr<MatMul1D> mat1_data;
    mat1_data.reset(buildThinStep2Matrix(ea,
                                         sig_sequence[dim],
                                         local_reps[dim],
                                         dim,
                                         m / mvec[dim],
                                         invert,
                                         /*inflate=*/true));
    matvec[dim].reset(new MatMul1DExec(*mat1_data, minimal));

    // Tower benchmark for the 97-component (D=96)
    if (local_reps[dim].length() == 96 && !invert) {
      std::cerr << "\n=== Tower Decomposition Analysis (dim=" << dim
                << ", D=96, mvec[dim]=" << mvec[dim] << ") ===\n";
      towerBenchmark96(ea, local_reps[dim], dim, m / mvec[dim]);
      std::cerr << "=== End Tower Analysis ===\n\n";
    }
  }

  for (long dim = nfactors - 2; dim >= 0; --dim) {
    std::unique_ptr<MatMul1D> mat_data;

    // Tower optimization: use butterfly factorization for D=96 dimension
    bool use_tower = (local_reps[dim].length() == 96 || local_reps[dim].length() == 256);

    if (use_tower && local_reps[dim].length() == 256) {
      // D=256: use standard ThinStep2Matrix (correct Vandermonde, BSGS optimized)
      std::cerr << "\n=== D=256: standard Vandermonde (correct) ===\n";
      mat_data.reset(buildThinStep2Matrix(ea, sig_sequence[dim],
                                          local_reps[dim], dim,
                                          m/mvec[dim], invert));
      matvec[dim].reset(new MatMul1DExec(*mat_data, minimal));
      std::cerr << "=== End D=256 ===\n\n";
    } else if (use_tower && local_reps[dim].length() == 96) {
      // D=96: Rader FACTORED - 4 multiplicative stages applied independently
      // V = P_out_DC * INTT * diag(B) * conj_NTT
      std::cerr << "\n=== TOWER D=96: Rader 4 INDEPENDENT stages ===\n";
      long D = 96;

      std::vector<NTL::Mat<NTL::zz_pX>> stages;
      buildRaderFactored96<PA_zz_p>(ea, local_reps[dim], dim,
                                     m / mvec[dim], invert, stages);

      long num_stages = stages.size();  // 4
      long old_len = matvec.length();
      matvec.SetLength(old_len + num_stages - 1);
      for (long i = old_len + num_stages - 2; i > dim + num_stages - 1; i--)
        matvec[i] = std::move(matvec[i - num_stages + 1]);

      for (long stage = 0; stage < num_stages; stage++) {
        long matvec_idx = dim + stage;
        std::unique_ptr<MatMul1D> stage_mat(
          new ThinStep2Tower96Stage<PA_zz_p>(ea, dim, stages[stage]));
        matvec[matvec_idx].reset(new MatMul1DExec(*stage_mat, minimal));
      }
      std::cerr << "  Installed " << num_stages << " independent Rader stages\n";
      std::cerr << "  (conj_NTT[96 diags] + diag(B)[1 diag] + INTT[96 diags] + P_out[96 diags])\n";
      std::cerr << "=== End Tower D=96 ===\n\n";
    } else {
      mat_data.reset(buildThinStep2Matrix(ea,
                                          sig_sequence[dim],
                                          local_reps[dim],
                                          dim,
                                          m / mvec[dim],
                                          invert));
      matvec[dim].reset(new MatMul1DExec(*mat_data, minimal));
    }
  }

  if (build_cache)
    upgrade();
}

void ThinEvalMap::upgrade()
{
  for (long i = 0; i < matvec.length(); i++)
    if (matvec[i])
      matvec[i]->upgrade();
}

// Applying the evaluation (or its inverse) map to a ciphertext
void ThinEvalMap::apply(Ctxt& ctxt) const
{
  if (!invert) { // forward direction
    for (long i = matvec.length() - 1; i >= 0; i--)
      if (matvec[i])
        matvec[i]->mul(ctxt);
  } else { // inverse transformation
    for (long i = 0; i < matvec.length(); i++)
      matvec[i]->mul(ctxt);
    traceMap(ctxt);
  }
}

// The callback interface for the matrix-multiplication routines.

//! \cond FALSE (make doxygen ignore these classes)
template <typename type>
class ThinStep2Matrix : public MatMul1D_derived<type>
{
  PA_INJECT(type)

  const EncryptedArray& base_ea;
  std::shared_ptr<CubeSignature> sig;
  long dim;
  NTL::Mat<RX> A;

public:
  // constructor
  ThinStep2Matrix(const EncryptedArray& _ea,
                  std::shared_ptr<CubeSignature> _sig,
                  const NTL::Vec<long>& reps,
                  long _dim,
                  long cofactor,
                  bool invert,
                  bool inflate) :
      base_ea(_ea), sig(_sig), dim(_dim)
  {
    long sz = sig->getDim(dim);
    assertEq(sz,
             reps.length(),
             "Invalid argument: sig and reps have inconsistent "
             "dimension");

    const EncryptedArrayDerived<type>& ea = _ea.getDerived(type());
    RBak bak;
    bak.save();
    _ea.getAlMod().restoreContext();
    const RX& G = ea.getG();
    long d = deg(G);

    NTL::Vec<RX> points(NTL::INIT_SIZE, sz);
    for (long j = 0; j < sz; j++) {
      points[j] = RX(reps[j] * cofactor, 1) % G;
      if (inflate)
        points[j] = NTL::PowerMod(points[j], d, G);
    }

    A.SetDims(sz, sz);
    for (long j = 0; j < sz; j++)
      A[0][j] = 1;

    for (long i = 1; i < sz; i++)
      for (long j = 0; j < sz; j++)
        A[i][j] = (A[i - 1][j] * points[j]) % G;

    if (invert) {
      REBak ebak;
      ebak.save();
      ea.restoreContextForG();

      mat_RE A1, A2;
      conv(A1, A);

      long p = _ea.getAlMod().getZMStar().getP();
      long r = _ea.getAlMod().getR();

      ppInvert(A2, A1, p, r);
      conv(A, A2);
    }
  }

  bool get(RX& out, long i, long j, UNUSED long k) const override
  {
    out = A[i][j];
    return false;
  }

  const EncryptedArray& getEA() const override { return base_ea; }
  bool multipleTransforms() const override { return false; }
  long getDim() const override { return dim; }
};

// ============================================================
// Tower 2-Stage Decomposition for D=96 Step2Matrix
// Stage 1: 6 diagonals at offsets {0,16,32,48,64,80}
// Stage 2: computed as S * Stage1^{-1}
// ============================================================

// (ThinStep2Tower96Stage class defined earlier, before constructor)

// Helper: build butterfly stage matrices with correct omega_3
template <typename type>
static void buildButterflyStages96(const EncryptedArray& _ea,
                                    const NTL::Vec<long>& reps,
                                    long dim,
                                    long cofactor,
                                    bool invert,
                                    std::vector<NTL::Mat<typename type::RX>>& stages_out)
{
  PA_INJECT(type)
  (void)tag; (void)dim;
  const long D = 96;

  const EncryptedArrayDerived<type>& ea = _ea.getDerived(type());
  RBak bak; bak.save();
  _ea.getAlMod().restoreContext();
  const RX& G = ea.getG();
  long p_val = _ea.getAlMod().getZMStar().getP();

  // Set up extension field
  REBak ebak; ebak.save();
  ea.restoreContextForG();

  // Find omega_32 (scalar, exists in Z_p)
  long omega32 = 0;
  for (long g = 2; g < 100; g++) {
    long val = NTL::PowerMod(g, (p_val - 1) / 32, p_val);
    if (NTL::PowerMod(val, 16, p_val) != 1) { omega32 = val; break; }
  }

  // Find omega_3 in F_{p^6} using RE (extension field element)
  // omega_3 = gen^{(p^6-1)/3}
  NTL::ZZ p6 = NTL::power(NTL::ZZ(p_val), 6);
  NTL::ZZ exp3 = (p6 - 1) / 3;

  RE omega3_E;
  for (long trial = 2; trial < 200; trial++) {
    RX gen_poly;
    NTL::SetCoeff(gen_poly, 0, trial);
    NTL::SetCoeff(gen_poly, 1, 1);
    RE gen_E; conv(gen_E, gen_poly);
    omega3_E = NTL::power(gen_E, exp3);
    if (!IsOne(omega3_E)) break;
  }
  RX omega3_poly = rep(omega3_E);

  std::cerr << "  omega_32=" << omega32 << ", omega_3 found (deg="
            << NTL::deg(omega3_poly) << ")" << std::endl;

  // Build 6 butterfly stages
  long radices[] = {3, 2, 2, 2, 2, 2};
  long strides_arr[] = {32, 16, 8, 4, 2, 1};
  stages_out.resize(6);

  for (long stage = 0; stage < 6; stage++) {
    long r = radices[stage];
    long s = strides_arr[stage];
    long gs = r * s;
    long num_groups = D / gs;

    stages_out[stage].SetDims(D, D);

    for (long grp = 0; grp < num_groups; grp++) {
      long base = grp * gs;
      for (long a = 0; a < r; a++) {
        for (long j = 0; j < s; j++) {
          for (long b = 0; b < r; b++) {
            long row = base + a * s + j;
            long col = base + b * s + j;

            if (r == 2) {
              long exp1 = (long)(((long long)a * b * (D / r)) % D);
              long exp2 = (long)(((long long)j * b * (D / gs)) % D);
              long val1 = NTL::PowerMod(omega32, exp1 % 32, p_val);
              long val2 = NTL::PowerMod(omega32, exp2 % 32, p_val);
              long val = (long)(((long long)val1 * val2) % p_val);
              stages_out[stage][row][col] = RX(val, 0);
            } else {
              // Radix-3: omega_3^{a*b} * omega_32^{j*b*(D/gs)}
              long ab_mod3 = (a * b) % 3;
              RE omega3_pow = NTL::power(omega3_E, ab_mod3);
              RX entry = rep(omega3_pow);
              long scalar_exp = (long)(((long long)j * b * (D / gs)) % D);
              long scalar_tw = NTL::PowerMod(omega32, scalar_exp % 32, p_val);
              // Multiply polynomial by scalar
              entry *= scalar_tw;
              entry %= G;
              stages_out[stage][row][col] = entry;
            }
          }
        }
      }
    }
  }

  // If inverse, rebuild stages with inverse twiddles (omega^{-1}) and 1/r scaling
  // instead of expensive matrix inversion
  if (invert) {
    stages_out.clear();
    stages_out.resize(6);

    // omega_32^{-1}
    long omega32_inv = NTL::InvMod(omega32, p_val);
    // omega_3^{-1} = omega_3^2 (since omega_3^3 = 1)
    RE omega3_inv_E = omega3_E * omega3_E;

    // Build inverse stages in REVERSE order
    for (long stage = 0; stage < 6; stage++) {
      long r = radices[stage];
      long s = strides_arr[stage];
      long gs = r * s;
      long num_groups = D / gs;
      long r_inv = NTL::InvMod(r, p_val);

      stages_out[stage].SetDims(D, D);

      for (long grp = 0; grp < num_groups; grp++) {
        long base = grp * gs;
        for (long a = 0; a < r; a++) {
          for (long j = 0; j < s; j++) {
            for (long b = 0; b < r; b++) {
              long row = base + a * s + j;
              long col = base + b * s + j;

              if (r == 2) {
                // Inverse: (1/r) * omega_32^{-a*b*(D/r)} * omega_32^{-j*a*(D/gs)}
                long exp1 = (long)(((long long)a * b * (D / r)) % D);
                long exp2 = (long)(((long long)j * a * (D / gs)) % D);
                long val1 = NTL::PowerMod(omega32_inv, exp1 % 32, p_val);
                long val2 = NTL::PowerMod(omega32_inv, exp2 % 32, p_val);
                long val = (long)((((long long)val1 * val2) % p_val * r_inv) % p_val);
                stages_out[stage][row][col] = RX(val, 0);
              } else {
                // Inverse radix-3: (1/r) * omega_3^{-a*b} * omega_32^{-j*a*(D/gs)}
                long ab_mod3 = (a * b) % 3;
                // omega_3^{-ab} = omega_3^{3-ab mod 3} = omega_3_inv^{ab}
                RE omega3_inv_pow = NTL::power(omega3_inv_E, ab_mod3);
                long scalar_exp = (long)(((long long)j * a * (D / gs)) % D);
                long scalar_tw = NTL::PowerMod(omega32_inv, scalar_exp % 32, p_val);
                scalar_tw = (long)(((long long)scalar_tw * r_inv) % p_val);
                RX entry = rep(omega3_inv_pow);
                entry *= scalar_tw;
                entry %= G;
                stages_out[stage][row][col] = entry;
              }
            }
          }
        }
      }
    }
    // Reverse stage order for inverse
    std::reverse(stages_out.begin(), stages_out.end());
  }
}

// Build the Rader-verified Vandermonde matrix for D=96
// Uses omega_96 in F_{p^d} computed via HElib's actual slot polynomial G
template <typename type>
static void buildRaderVandermonde96(const EncryptedArray& _ea,
                                    const NTL::Vec<long>& reps,
                                    long dim, long cofactor, bool invert,
                                    NTL::Mat<typename type::RX>& V_out)
{
  PA_INJECT(type)
  (void)tag; (void)dim;
  const long D = 96;
  const long ell = 97;
  const long g_ell = 5;

  const EncryptedArrayDerived<type>& ea = _ea.getDerived(type());
  RBak bak; bak.save();
  _ea.getAlMod().restoreContext();
  const RX& G = ea.getG();
  long p_val = _ea.getAlMod().getZMStar().getP();

  REBak ebak; ebak.save();
  ea.restoreContextForG();

  // omega_96 in F_{p^d}
  NTL::ZZ pd = NTL::power(NTL::ZZ(p_val), deg(G));
  NTL::ZZ expD = (pd - 1) / D;
  RE omega_D;
  for (long trial = 2; trial < 200; trial++) {
    RX gen_poly; SetCoeff(gen_poly, 0, trial); SetCoeff(gen_poly, 1, 1);
    RE gen; conv(gen, gen_poly);
    omega_D = power(gen, expD);
    if (!IsOne(omega_D) && IsOne(power(omega_D, D))) {
      bool good = true;
      if (IsOne(power(omega_D, D/2))) good = false;
      if (IsOne(power(omega_D, D/3))) good = false;
      if (good) break;
    }
  }

  // alpha and evaluation points
  RX X_poly; SetCoeff(X_poly, 1, 1);
  RE X_elem; conv(X_elem, X_poly);
  RE alpha = power(X_elem, cofactor);

  NTL::Vec<RE> points;
  points.SetLength(D);
  for (long j = 0; j < D; j++) {
    long exp = 1;
    for (long i = 0; i < j; i++) exp = (exp * g_ell) % ell;
    points[j] = power(alpha, exp);
  }

  // Rader kernel: B[k] = sum_t point[t] * omega_D^{tk}
  NTL::Vec<RE> B;
  B.SetLength(D);
  for (long k = 0; k < D; k++) {
    clear(B[k]);
    for (long t = 0; t < D; t++) {
      long exp = (long)(((long long)t * k) % D);
      B[k] += points[t] * power(omega_D, exp);
    }
  }

  // dlog table
  std::vector<long> dlog(ell, -1);
  long g_pow = 1;
  for (long s = 0; s < D; s++) { dlog[g_pow] = s; g_pow = (g_pow * g_ell) % ell; }

  // Build V via Rader formula:
  // V[0][j] = 1
  // V[i][j] = (1/D) * sum_k B[k] * omega_D^{-(dlog[i]+j)*k}  for i=1,...,D-1
  RE omega_inv = power(omega_D, D - 1);
  long D_inv_scalar = NTL::InvMod(D, p_val);
  RX D_inv_poly; SetCoeff(D_inv_poly, 0, D_inv_scalar);
  RE D_inv_E; conv(D_inv_E, D_inv_poly);

  V_out.SetDims(D, D);
  for (long j = 0; j < D; j++) NTL::set(V_out[0][j]);
  for (long i = 1; i < D; i++) {
    long s = dlog[i];
    for (long j = 0; j < D; j++) {
      RE val; clear(val);
      for (long k = 0; k < D; k++) {
        long exp = (long)((long long)(s + j) * k % D);
        val += B[k] * power(omega_inv, exp);
      }
      val *= D_inv_E;
      V_out[i][j] = rep(val);
    }
  }

  // If inverse: invert V
  if (invert) {
    mat_RE V_E; conv(V_E, V_out);
    mat_RE V_inv;
    ppInvert(V_inv, V_E, p_val, _ea.getAlMod().getR());
    conv(V_out, V_inv);
  }

  std::cerr << "  [Rader] Built " << (invert ? "inverse " : "") << "Vandermonde via Rader (D="
            << D << ", omega_96 order verified)" << std::endl;
}

// Build Rader factored stages: returns 4 matrices whose PRODUCT = V
// stages[0] = conj_NTT, stages[1] = diag(B), stages[2] = INTT, stages[3] = P_out_DC
// For inverse: returns inverse factors in reverse order
template <typename type>
static void buildRaderFactored96(const EncryptedArray& _ea,
                                  const NTL::Vec<long>& reps,
                                  long dim, long cofactor, bool invert,
                                  std::vector<NTL::Mat<typename type::RX>>& stages_out)
{
  PA_INJECT(type)
  (void)tag; (void)dim;
  const long D = 96;
  const long ell = 97;
  const long g_ell = 5;

  const EncryptedArrayDerived<type>& ea = _ea.getDerived(type());
  RBak bak; bak.save();
  _ea.getAlMod().restoreContext();
  const RX& G = ea.getG();
  long p_val = _ea.getAlMod().getZMStar().getP();

  REBak ebak; ebak.save();
  ea.restoreContextForG();

  // Compute omega_96
  NTL::ZZ pd = NTL::power(NTL::ZZ(p_val), deg(G));
  NTL::ZZ expD = (pd - 1) / D;
  RE omega_D;
  for (long trial = 2; trial < 200; trial++) {
    RX gen_poly; SetCoeff(gen_poly, 0, trial); SetCoeff(gen_poly, 1, 1);
    RE gen; conv(gen, gen_poly);
    omega_D = power(gen, expD);
    if (!IsOne(omega_D) && IsOne(power(omega_D, D))) {
      bool good = true;
      if (IsOne(power(omega_D, D/2))) good = false;
      if (IsOne(power(omega_D, D/3))) good = false;
      if (good) break;
    }
  }
  RE omega_inv = power(omega_D, D - 1);

  // alpha and points
  RX X_poly; SetCoeff(X_poly, 1, 1);
  RE X_elem; conv(X_elem, X_poly);
  RE alpha = power(X_elem, cofactor);

  NTL::Vec<RE> points;
  points.SetLength(D);
  for (long j = 0; j < D; j++) {
    long exp = 1;
    for (long i = 0; i < j; i++) exp = (exp * g_ell) % ell;
    points[j] = power(alpha, exp);
  }

  // Kernel B[k] = sum_t point[t] * omega^{tk}
  NTL::Vec<RE> B;
  B.SetLength(D);
  for (long k = 0; k < D; k++) {
    clear(B[k]);
    for (long t = 0; t < D; t++) {
      long exp = (long)(((long long)t * k) % D);
      B[k] += points[t] * power(omega_D, exp);
    }
  }

  // dlog table
  std::vector<long> dlog(ell, -1);
  long g_pow = 1;
  for (long s = 0; s < D; s++) { dlog[g_pow] = s; g_pow = (g_pow * g_ell) % ell; }

  long D_inv_scalar = NTL::InvMod(D, p_val);
  RX D_inv_poly; SetCoeff(D_inv_poly, 0, D_inv_scalar);
  RE D_inv_E; conv(D_inv_E, D_inv_poly);

  stages_out.resize(4);

  // Stage 0: conj_NTT[k][j] = omega^{-jk}
  stages_out[0].SetDims(D, D);
  for (long k = 0; k < D; k++)
    for (long j = 0; j < D; j++) {
      long exp = (long)(((long long)j * k) % D);
      stages_out[0][k][j] = rep(power(omega_inv, exp));
    }

  // Stage 1: diag(B)[k][k] = B[k]
  stages_out[1].SetDims(D, D);
  for (long k = 0; k < D; k++)
    stages_out[1][k][k] = rep(B[k]);

  // Stage 2: INTT[s][k] = (1/D) * omega^{-sk}
  stages_out[2].SetDims(D, D);
  for (long s = 0; s < D; s++)
    for (long k = 0; k < D; k++) {
      long exp = (long)(((long long)s * k) % D);
      RE val = power(omega_inv, exp) * D_inv_E;
      stages_out[2][s][k] = rep(val);
    }

  // Stage 3: P_out_DC[i][s] = 1 if i = g^s mod ell (for i>=1), row 0 = all ones
  stages_out[3].SetDims(D, D);
  for (long j = 0; j < D; j++) NTL::set(stages_out[3][0][j]);  // row 0: DC
  for (long i = 1; i < D; i++) {
    if (i < ell && dlog[i] >= 0)
      NTL::set(stages_out[3][i][dlog[i]]);  // permutation
  }

  // For inverse: invert each stage and reverse order
  if (invert) {
    // Stage 0 inv: (1/D) * NTT (omega, not omega_inv)
    for (long k = 0; k < D; k++)
      for (long j = 0; j < D; j++) {
        long exp = (long)(((long long)j * k) % D);
        RE val = power(omega_D, exp) * D_inv_E;
        stages_out[0][k][j] = rep(val);
      }
    // Stage 1 inv: diag(1/B[k])
    for (long k = 0; k < D; k++) {
      RE B_inv = power(B[k], pd - NTL::ZZ(2));  // B[k]^{-1} = B[k]^{p^d-2}
      stages_out[1][k][k] = rep(B_inv);
    }
    // Stage 2 inv: NTT (omega, not omega_inv), no 1/D scaling
    for (long s = 0; s < D; s++)
      for (long k = 0; k < D; k++) {
        long exp = (long)(((long long)s * k) % D);
        stages_out[2][s][k] = rep(power(omega_D, exp));
      }
    // Stage 3 inv: P_out^{-1} = P_out^T (transpose of permutation)
    NTL::Mat<RX> P_inv;
    P_inv.SetDims(D, D);
    for (long j = 0; j < D; j++) NTL::set(P_inv[0][j]);  // row 0 stays
    for (long i = 1; i < D; i++) {
      if (i < ell && dlog[i] >= 0)
        NTL::set(P_inv[dlog[i]][i]);  // transpose
    }
    stages_out[3] = P_inv;

    // Reverse order for inverse: V^{-1} = stage0_inv * stage1_inv * stage2_inv * stage3_inv
    std::reverse(stages_out.begin(), stages_out.end());
  }

  std::cerr << "  [Rader] Built 4 factored stages"
            << (invert ? " (inverse)" : " (forward)") << std::endl;
}

template <typename type>
static void buildTowerStages96(const EncryptedArray& _ea,
                               const NTL::Vec<long>& reps,
                               long dim,
                               long cofactor,
                               NTL::Mat<typename type::RX>& M1_out,
                               NTL::Mat<typename type::RX>& M2_out,
                               long& m2_ndiags)
{
  PA_INJECT(type)
  (void)tag; // suppress unused warning from PA_INJECT
  (void)dim;
  const long D = 96;

  const EncryptedArrayDerived<type>& ea = _ea.getDerived(type());
  RBak bak; bak.save();
  _ea.getAlMod().restoreContext();
  const RX& G = ea.getG();

  // Build evaluation points
  NTL::Vec<RX> points(NTL::INIT_SIZE, D);
  for (long j = 0; j < D; j++)
    points[j] = RX(reps[j] * cofactor, 1) % G;

  // Build full Step2Matrix S (Vandermonde)
  NTL::Mat<RX> S;
  S.SetDims(D, D);
  for (long j = 0; j < D; j++)
    set(S[0][j]);
  for (long i = 1; i < D; i++)
    for (long j = 0; j < D; j++)
      MulMod(S[i][j], S[i-1][j], points[j], G);

  // Build Stage 1: block-diagonal 6-point Vandermonde within each orbit
  // Orbit b = {b, b+16, b+32, b+48, b+64, b+80} for b=0,...,15
  // M1[b+16*r][b+16*c] = points[b+16*c]^r
  M1_out.SetDims(D, D);
  for (long b = 0; b < 16; b++) {
    for (long r = 0; r < 6; r++) {
      for (long c = 0; c < 6; c++) {
        long row = b + 16 * r;
        long col = b + 16 * c;
        PowerMod(M1_out[row][col], points[col], r, G);
      }
    }
  }

  // Invert M1 block-by-block using RE (extension field) arithmetic
  NTL::Mat<RX> M1_inv;
  M1_inv.SetDims(D, D);

  {
    REBak ebak; ebak.save();
    ea.restoreContextForG();

    for (long b = 0; b < 16; b++) {
      mat_RE block_E;
      block_E.SetDims(6, 6);
      for (long r = 0; r < 6; r++)
        for (long c = 0; c < 6; c++)
          conv(block_E[r][c], M1_out[b + 16*r][b + 16*c]);

      mat_RE inv_E;
      NTL::inv(inv_E, block_E);

      for (long r = 0; r < 6; r++)
        for (long c = 0; c < 6; c++)
          conv(M1_inv[b + 16*r][b + 16*c], inv_E[r][c]);
    }
  }

  // Compute M2 = S * M1_inv (mod G)
  M2_out.SetDims(D, D);
  RX tmp;
  for (long i = 0; i < D; i++) {
    for (long j = 0; j < D; j++) {
      clear(M2_out[i][j]);
      for (long l = 0; l < D; l++) {
        if (IsZero(S[i][l]) || IsZero(M1_inv[l][j])) continue;
        MulMod(tmp, S[i][l], M1_inv[l][j], G);
        M2_out[i][j] += tmp;
      }
      M2_out[i][j] %= G;
    }
  }

  // Count M2 non-zero diagonals
  std::vector<bool> seen(D, false);
  for (long i = 0; i < D; i++)
    for (long j = 0; j < D; j++)
      if (!IsZero(M2_out[i][j]))
        seen[((j - i) % D + D) % D] = true;
  m2_ndiags = 0;
  for (long i = 0; i < D; i++)
    if (seen[i]) m2_ndiags++;
}

// Public entry point for tower benchmark
static void towerBenchmark96(const EncryptedArray& ea,
                      const NTL::Vec<long>& reps,
                      long dim,
                      long cofactor)
{
  using namespace NTL;
  const long D = 96;

  std::cerr << "  [Tower] Building stages for D=" << D << "...\n";

  Mat<zz_pX> M1, M2;
  long m2_ndiags = 0;
  buildTowerStages96<PA_zz_p>(ea, reps, dim, cofactor, M1, M2, m2_ndiags);

  std::cerr << "  [Tower] Stage2 non-zero diagonals: " << m2_ndiags << "\n";

  // Count Stage1 diagonals
  std::vector<bool> seen1(D, false);
  for (long i = 0; i < D; i++)
    for (long j = 0; j < D; j++)
      if (!IsZero(M1[i][j]))
        seen1[((j-i)%D+D)%D] = true;
  long m1_ndiags = 0;
  for (long i = 0; i < D; i++) if (seen1[i]) m1_ndiags++;
  std::cerr << "  [Tower] Stage1 non-zero diagonals: " << m1_ndiags << "\n";

  if (m2_ndiags <= 16) {
    std::cerr << "  [Tower] *** SUCCESS: Tower factorization is valid! ***\n";
    std::cerr << "  [Tower] Cost: " << m1_ndiags << " + " << m2_ndiags
              << " = " << (m1_ndiags + m2_ndiags) << " diagonals\n";
    std::cerr << "  [Tower] vs full: " << D << " diagonals\n";
  } else {
    std::cerr << "  [Tower] Stage2 has " << m2_ndiags
              << " diags (expected <=16). Factorization NOT sparse.\n";
  }
}

// ============================================================

static MatMul1D* buildThinStep2Matrix(const EncryptedArray& ea,
                                      std::shared_ptr<CubeSignature> sig,
                                      const NTL::Vec<long>& reps,
                                      long dim,
                                      long cofactor,
                                      bool invert,
                                      bool inflate)
{
  switch (ea.getTag()) {
  case PA_GF2_tag:
    return new ThinStep2Matrix<PA_GF2>(ea,
                                       sig,
                                       reps,
                                       dim,
                                       cofactor,
                                       invert,
                                       inflate);

  case PA_zz_p_tag:
    return new ThinStep2Matrix<PA_zz_p>(ea,
                                        sig,
                                        reps,
                                        dim,
                                        cofactor,
                                        invert,
                                        inflate);

  default:
    return 0;
  }
}

template <typename type>
class ThinStep1Matrix : public MatMul1D_derived<type>
{
  PA_INJECT(type)

  const EncryptedArray& base_ea;
  std::shared_ptr<CubeSignature> sig;
  long dim;
  NTL::Mat<RX> A_deflated;

public:
  // constructor
  ThinStep1Matrix(const EncryptedArray& _ea,
                  std::shared_ptr<CubeSignature> _sig,
                  const NTL::Vec<long>& reps,
                  long _dim,
                  long cofactor) :
      base_ea(_ea), sig(_sig), dim(_dim)
  {
    const EncryptedArrayDerived<type>& ea = _ea.getDerived(type());
    RBak bak;
    bak.save();
    _ea.getAlMod().restoreContext();
    const RXModulus G(ea.getG());
    long d = deg(G);

    long p = _ea.getAlMod().getZMStar().getP();
    long r = _ea.getAlMod().getR();

    long sz = sig->getDim(dim);
    assertEq(sz,
             reps.length(),
             "Invalid argument: sig and reps have inconsistent "
             "dimension");
    assertEq(dim,
             sig->getNumDims() - 1,
             "Invalid argument: dim must be one less than "
             "sig->getNumDims()");
    assertEq(sig->getSize(), ea.size(), "sig and ea do not have matching size");

    // so sz == phi(m_last)/d, where d = deg(G) = order of p mod m

    NTL::Vec<RX> points(NTL::INIT_SIZE, sz);
    for (long j = 0; j < sz; j++)
      points[j] = RX(reps[j] * cofactor, 1) % G;

    NTL::Mat<RX> AA(NTL::INIT_SIZE, sz * d, sz);
    for (long j = 0; j < sz; j++)
      AA[0][j] = 1;

    for (long i = 1; i < sz * d; i++)
      for (long j = 0; j < sz; j++)
        AA[i][j] = (AA[i - 1][j] * points[j]) % G;

    NTL::Mat<mat_R> A;
    A.SetDims(sz, sz);
    for (long i = 0; i < sz; i++)
      for (long j = 0; j < sz; j++) {
        A[i][j].SetDims(d, d);
        for (long k = 0; k < d; k++)
          VectorCopy(A[i][j][k], AA[i * d + k][j], d);
      }

    // if (invert) {
    // NOTE: this version is only used for the inverse matrix
    mat_R A1, A2;
    A1.SetDims(sz * d, sz * d);
    for (long i = 0; i < sz * d; i++)
      for (long j = 0; j < sz * d; j++)
        A1[i][j] = A[i / d][j / d][i % d][j % d];

    ppInvert(A2, A1, p, r);

    for (long i = 0; i < sz * d; i++)
      for (long j = 0; j < sz * d; j++)
        A[i / d][j / d][i % d][j % d] = A2[i][j];
    // }

    A_deflated.SetDims(sz, sz);
    vec_R v, w;
    v.SetLength(d);
    w.SetLength(d);

    RX h; // set h = X^p mod G
    PowerXMod(h, p, G);

    NTL::Vec<R> trace_vec;
    trace_vec.SetLength(2 * d - 1);
    // set trace_vec[i] = trace(X^i mod G)
    for (long i = 0; i < 2 * d - 1; i++) {
      RX trace_val;
      TraceMap(trace_val, (RX(i, 1) % G), d, G, h);
      assertTrue(deg(trace_val) <= 0, "trace_val is positive");
      trace_vec[i] = ConstTerm(trace_val);
    }

    NTL::Mat<R> trace_mat;
    trace_mat.SetDims(d, d);
    // set trace_mat[i][j] = trace(X^{i+j} mod G)
    for (long i = 0; i < d; i++)
      for (long j = 0; j < d; j++)
        trace_mat[i][j] = trace_vec[i + j];

    NTL::Mat<R> trace_mat_inv;
    RelaxedInv(trace_mat_inv, trace_mat);

    for (long i = 0; i < sz; i++)
      for (long j = 0; j < sz; j++) {
        for (long i1 = 0; i1 < d; i1++)
          v[i1] = A[i][j][i1][0];
        mul(w, v, trace_mat_inv);
        conv(A_deflated[i][j], w);
      }
  } // constructor

  bool get(RX& out, long i, long j, UNUSED long k) const override
  {
    out = A_deflated[i][j];
    return false;
  }

  const EncryptedArray& getEA() const override { return base_ea; }
  bool multipleTransforms() const override { return false; }
  long getDim() const override { return dim; }
};

static MatMul1D* buildThinStep1Matrix(const EncryptedArray& ea,
                                      std::shared_ptr<CubeSignature> sig,
                                      const NTL::Vec<long>& reps,
                                      long dim,
                                      long cofactor)
{
  switch (ea.getTag()) {
  case PA_GF2_tag:
    return new ThinStep1Matrix<PA_GF2>(ea, sig, reps, dim, cofactor);

  case PA_zz_p_tag:
    return new ThinStep1Matrix<PA_zz_p>(ea, sig, reps, dim, cofactor);

  default:
    return 0;
  }
}
//! \endcond

} // namespace helib
