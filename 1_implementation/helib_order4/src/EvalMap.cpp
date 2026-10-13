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

#include <cstdlib>
#include <fstream>
#include <sstream>

// needed to get NTL's TraceMap functions...needed for ThinEvalMap
#include <NTL/lzz_pXFactoring.h>
#include <NTL/GF2XFactoring.h>

#include "tensor_T16_coeffs.h"

namespace helib {

// Forward declarations
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
static MatMulExecBase* buildThinStep2Rader97Exec(const EncryptedArray& ea,
                                                 std::shared_ptr<CubeSignature> sig,
                                                 const NTL::Vec<long>& reps,
                                                 long dim,
                                                 long cofactor,
                                                 bool invert,
                                                 bool inflate);
static MatMulExecBase* buildThinStep2Radix4Exec(const EncryptedArray& ea,
                                                std::shared_ptr<CubeSignature> sig,
                                                const NTL::Vec<long>& reps,
                                                long dim,
                                                long cofactor,
                                                bool invert,
                                                bool inflate);
static MatMulExecBase* buildThinStep2Radix4Po2Exec(const EncryptedArray& ea,
                                                   std::shared_ptr<CubeSignature> sig,
                                                   const NTL::Vec<long>& reps,
                                                   long dim,
                                                   long cofactor,
                                                   bool invert,
                                                   bool inflate);
static MatMulExecBase* buildThinStep2TensorExec(const EncryptedArray& ea,
                                                std::shared_ptr<CubeSignature> sig,
                                                const NTL::Vec<long>& reps,
                                                long dim,
                                                long cofactor,
                                                bool invert,
                                                bool inflate);
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

    const EncryptedArrayDerived<type>& ea UNUSED = _ea.getDerived(type());
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

    const char* dump_flag = std::getenv("HELIB_DUMP_STEP2_INTERNAL");
    if (dump_flag) {
      std::string path = std::string("step2_matrix_dim") + std::to_string(dim) + "_sz" + std::to_string(sz) + ".txt";
      std::ofstream ofs(path);
      if (ofs.is_open()) {
        ofs << sz << "\n";
        for (long i = 0; i < sz; i++) {
          for (long j = 0; j < sz; j++) {
            std::ostringstream ss;
            ss << A[i][j];
            ofs << ss.str();
            if (j < sz - 1) ofs << "|";
          }
          ofs << "\n";
        }
        ofs.close();
        std::cout << "Step2Matrix dumped: " << path << " (" << sz << "x" << sz << ")\n";
      }
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

// Fat Step2Matrix ends here

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
    const EncryptedArrayDerived<type>& ea UNUSED = _ea.getDerived(type());
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
    if (MatMulExecBase* exec =
            buildThinStep2Radix4Po2Exec(ea,
                                        sig_sequence[dim],
                                        local_reps[dim],
                                        dim,
                                        m / mvec[dim],
                                        invert,
                                        /*inflate=*/true)) {
      matvec[dim].reset(exec);
    } else if (MatMulExecBase* exec =
            buildThinStep2Radix4Exec(ea,
                                     sig_sequence[dim],
                                     local_reps[dim],
                                     dim,
                                     m / mvec[dim],
                                     invert,
                                     /*inflate=*/true)) {
      matvec[dim].reset(exec);
    } else if (MatMulExecBase* exec =
            buildThinStep2Rader97Exec(ea,
                                      sig_sequence[dim],
                                      local_reps[dim],
                                      dim,
                                      m / mvec[dim],
                                      invert,
                                      /*inflate=*/true)) {
      matvec[dim].reset(exec);
    } else {
      std::unique_ptr<MatMul1D> mat1_data;
      mat1_data.reset(buildThinStep2Matrix(ea,
                                           sig_sequence[dim],
                                           local_reps[dim],
                                           dim,
                                           m / mvec[dim],
                                           invert,
                                           /*inflate=*/true));
      matvec[dim].reset(new MatMul1DExec(*mat1_data, minimal));
    }
  }

  for (long dim = nfactors - 2; dim >= 0; --dim) {
    if (MatMulExecBase* exec =
            buildThinStep2Radix4Po2Exec(ea,
                                        sig_sequence[dim],
                                        local_reps[dim],
                                        dim,
                                        m / mvec[dim],
                                        invert,
                                        /*inflate=*/false)) {
      matvec[dim].reset(exec);
    } else if (MatMulExecBase* exec =
            buildThinStep2Radix4Exec(ea,
                                     sig_sequence[dim],
                                     local_reps[dim],
                                     dim,
                                     m / mvec[dim],
                                     invert,
                                     /*inflate=*/false)) {
      matvec[dim].reset(exec);
    } else if (MatMulExecBase* exec =
            buildThinStep2TensorExec(ea,
                                     sig_sequence[dim],
                                     local_reps[dim],
                                     dim,
                                     m / mvec[dim],
                                     invert,
                                     /*inflate=*/false)) {
      matvec[dim].reset(exec);
    } else if (MatMulExecBase* exec =
            buildThinStep2Rader97Exec(ea,
                                      sig_sequence[dim],
                                      local_reps[dim],
                                      dim,
                                      m / mvec[dim],
                                      invert,
                                      /*inflate=*/false)) {
      matvec[dim].reset(exec);
    } else {
      std::unique_ptr<MatMul1D> mat_data;

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

void ThinEvalMap::apply(Ctxt& ctxt0, Ctxt& ctxt1) const
{
  if (!invert) { // forward direction
    for (long i = matvec.length() - 1; i >= 0; i--) {
      if (matvec[i]) {
        matvec[i]->mul(ctxt0);
        matvec[i]->mul(ctxt1);
      }
    }
  } else { // inverse transformation
    for (long i = 0; i < matvec.length(); i++) {
      if (matvec[i]) {
        matvec[i]->mul(ctxt0);
        matvec[i]->mul(ctxt1);
      }
    }
    traceMap(ctxt0);
    traceMap(ctxt1);
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

    const EncryptedArrayDerived<type>& ea UNUSED = _ea.getDerived(type());
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

    const char* dump_flag = std::getenv("HELIB_DUMP_THINSTEP2");
    if (dump_flag) {
      std::string path = std::string("thinstep2_dim") + std::to_string(dim) +
                         "_sz" + std::to_string(sz) + "_d" + std::to_string(d) +
                         (inflate ? "_inflated" : "") +
                         (invert ? "_inv" : "") + ".txt";
      std::ofstream ofs(path);
      if (ofs.is_open()) {
        ofs << "sz " << sz << "\n";
        ofs << "d " << d << "\n";
        ofs << "inflate " << (inflate ? 1 : 0) << "\n";
        ofs << "invert " << (invert ? 1 : 0) << "\n";
        for (long i = 0; i < sz; i++) {
          for (long j = 0; j < sz; j++) {
            std::ostringstream ss;
            ss << A[i][j];
            ofs << ss.str();
            if (j < sz - 1) ofs << "|";
          }
          ofs << "\n";
        }
        ofs.close();
        std::cout << "ThinStep2Matrix dumped: " << path
                  << " (sz=" << sz << ", d=" << d << ", inflate=" << inflate << ")\n";
      }
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

static bool thinStep2Rader97Enabled()
{
  const char* flag = std::getenv("HELIB_THINSTEP2_RADER97");
  return flag != nullptr && std::atol(flag) != 0;
}

static bool thinStep2Radix4Enabled()
{
  const char* flag = std::getenv("HELIB_THINSTEP2_RADIX4");
  return flag != nullptr && std::atol(flag) != 0;
}

template <typename type>
class Dense1DMatrix : public MatMul1D_derived<type>
{
  PA_INJECT(type)

  const EncryptedArray& base_ea;
  long dim;
  NTL::Mat<RX> A;

public:
  Dense1DMatrix(const EncryptedArray& _ea, long _dim, const NTL::Mat<RX>& _A) :
      base_ea(_ea), dim(_dim), A(_A)
  {}

  bool get(RX& out, long i, long j, UNUSED long k) const override
  {
    out = A[i][j];
    return IsZero(out);
  }

  const EncryptedArray& getEA() const override { return base_ea; }
  bool multipleTransforms() const override { return false; }
  long getDim() const override { return dim; }
};

template <typename type>
static typename type::RX rxConst(long value)
{
  using R = typename type::R;
  using RX = typename type::RX;
  RX out;
  clear(out);
  SetCoeff(out, 0, R(value));
  return out;
}

template <typename type>
static typename type::RX rxConst(const typename type::R& value)
{
  using RX = typename type::RX;
  RX out;
  clear(out);
  SetCoeff(out, 0, value);
  return out;
}

template <typename type>
static typename type::RX rxMulMod(const typename type::RX& a,
                                  const typename type::RX& b,
                                  const typename type::RX& G)
{
  using RX = typename type::RX;
  RX out;
  MulMod(out, a, b, G);
  return out;
}

template <typename type>
static typename type::RX rxPowerMod(const typename type::RX& a,
                                    long e,
                                    const typename type::RX& G)
{
  return NTL::PowerMod(a, e, G);
}

template <typename type>
static std::vector<long> diagonalOffsets(const NTL::Mat<typename type::RX>& A)
{
  long D = A.NumRows();
  std::vector<bool> seen(D, false);
  for (long i : range(D))
    for (long j : range(D))
      if (!IsZero(A[i][j]))
        seen[mcMod(j - i, D)] = true;

  std::vector<long> offsets;
  for (long i : range(D))
    if (seen[i])
      offsets.push_back(i);
  return offsets;
}

template <typename type>
static void addStage(std::vector<std::unique_ptr<SparseMatMul1DExec>>& stages,
                     const EncryptedArray& ea,
                     long dim,
                     const NTL::Mat<typename type::RX>& A)
{
  Dense1DMatrix<type> mat(ea, dim, A);
  stages.emplace_back(new SparseMatMul1DExec(mat, diagonalOffsets<type>(A)));
}

template <typename type>
static NTL::Mat<typename type::RX> zeroMatrix(long D)
{
  using RX = typename type::RX;
  NTL::Mat<RX> A;
  A.SetDims(D, D);
  return A;
}

template <typename type>
static NTL::Mat<typename type::RX>
buildRaderInputCRTMatrix(const NTL::Vec<long>& reps)
{
  using RX = typename type::RX;
  const long D = 96;
  NTL::Mat<RX> A = zeroMatrix<type>(D);
  RX one = rxConst<type>(1);
  RX minus_one = rxConst<type>(-1);

  for (long t : range(D)) {
    long col = (t % 3) * 32 + (t % 32);
    A[0][col] += minus_one;

    long coeff_idx = reps[mcMod(-t, D)];
    if (coeff_idx < D)
      A[coeff_idx][col] += one;
  }

  return A;
}

template <typename type>
static NTL::Mat<typename type::RX>
buildCRTOutputMatrix()
{
  using RX = typename type::RX;
  const long D = 96;
  NTL::Mat<RX> A = zeroMatrix<type>(D);
  RX one = rxConst<type>(1);

  for (long t : range(D)) {
    long row = (t % 3) * 32 + (t % 32);
    A[row][t] = one;
  }

  return A;
}

template <typename type>
static NTL::Mat<typename type::RX>
buildDIF32StageMatrix(long n,
                      const typename type::RX& root,
                      const typename type::RX& G)
{
  using RX = typename type::RX;
  const long D = 96;
  long m = n / 2;
  NTL::Mat<RX> A = zeroMatrix<type>(D);
  RX omega = rxPowerMod<type>(root, m, G);

  for (long block : range(3)) {
    long block_base = block * 32;
    for (long seg = 0; seg < 32; seg += n) {
      long base = block_base + seg;
      for (long q : range(m)) {
        for (long s : range(2)) {
          long in = base + q + m * s;
          for (long k : range(2)) {
            long out = base + k * m + q;
            RX val =
                rxMulMod<type>(rxPowerMod<type>(root, q * k, G),
                               rxPowerMod<type>(omega, s * k, G),
                               G);
            A[in][out] += val;
          }
        }
      }
    }
  }

  return A;
}

template <typename type>
static NTL::Mat<typename type::RX>
buildDIF32Radix4StageMatrix(long n,
                            const typename type::RX& root,
                            const typename type::RX& G)
{
  using RX = typename type::RX;
  const long D = 96;
  long m = n / 4;
  NTL::Mat<RX> A = zeroMatrix<type>(D);
  RX omega4 = rxPowerMod<type>(root, m, G);

  for (long block : range(3)) {
    long block_base = block * 32;
    for (long seg = 0; seg < 32; seg += n) {
      long base = block_base + seg;
      for (long q : range(m)) {
        for (long s : range(4)) {
          long in = base + q + m * s;
          for (long k : range(4)) {
            long out = base + k * m + q;
            RX val =
                rxMulMod<type>(rxPowerMod<type>(root, q * k, G),
                               rxPowerMod<type>(omega4, s * k, G),
                               G);
            A[in][out] += val;
          }
        }
      }
    }
  }

  return A;
}

template <typename type>
static NTL::Mat<typename type::RX>
buildInverseDIF32Radix4StageMatrix(long n,
                                   const typename type::RX& root,
                                   const typename type::RX& G)
{
  using R = typename type::R;
  using RX = typename type::RX;
  const long D = 96;
  long m = n / 4;
  NTL::Mat<RX> A = zeroMatrix<type>(D);
  RX omega4 = rxPowerMod<type>(root, m, G);
  R inv4_R = inv(R(4));
  RX inv4 = rxConst<type>(inv4_R);

  for (long block : range(3)) {
    long block_base = block * 32;
    for (long seg = 0; seg < 32; seg += n) {
      long base = block_base + seg;
      for (long q : range(m)) {
        for (long s : range(4)) {
          long out = base + q + m * s;
          for (long k : range(4)) {
            long in = base + k * m + q;
            RX val = rxPowerMod<type>(root, mcMod(-q * k, n), G);
            val = rxMulMod<type>(val,
                                 rxPowerMod<type>(omega4, mcMod(-s * k, 4), G),
                                 G);
            val = rxMulMod<type>(val, inv4, G);
            A[in][out] += val;
          }
        }
      }
    }
  }

  return A;
}

template <typename type>
static NTL::Mat<typename type::RX>
buildInverseDIF32StageMatrix(long n,
                             const typename type::RX& root,
                             const typename type::RX& G)
{
  using R = typename type::R;
  using RX = typename type::RX;
  const long D = 96;
  long m = n / 2;
  NTL::Mat<RX> A = zeroMatrix<type>(D);
  RX omega = rxPowerMod<type>(root, m, G);
  R inv2_R = inv(R(2));
  RX inv2 = rxConst<type>(inv2_R);

  for (long block : range(3)) {
    long block_base = block * 32;
    for (long seg = 0; seg < 32; seg += n) {
      long base = block_base + seg;
      for (long q : range(m)) {
        for (long s : range(2)) {
          long out = base + q + m * s;
          for (long k : range(2)) {
            long in = base + k * m + q;
            RX val = rxPowerMod<type>(root, mcMod(-q * k, n), G);
            val = rxMulMod<type>(val,
                                 rxPowerMod<type>(omega, mcMod(-s * k, 2), G),
                                 G);
            val = rxMulMod<type>(val, inv2, G);
            A[in][out] += val;
          }
        }
      }
    }
  }

  return A;
}

template <typename type>
static void applyDIF32ToBlock(std::vector<typename type::RX>& vec,
                              const typename type::RX& root32,
                              const typename type::RX& G)
{
  using RX = typename type::RX;
  for (long n = 32; n >= 2; n /= 2) {
    long m = n / 2;
    RX root = rxPowerMod<type>(root32, 32 / n, G);
    RX omega = rxPowerMod<type>(root, m, G);
    std::vector<RX> next = vec;

    for (long seg = 0; seg < 32; seg += n) {
      for (long q : range(m)) {
        for (long k : range(2)) {
          RX acc;
          clear(acc);
          for (long s : range(2)) {
            RX val = rxMulMod<type>(vec[seg + q + m * s],
                                    rxPowerMod<type>(root, q * k, G),
                                    G);
            val = rxMulMod<type>(val,
                                 rxPowerMod<type>(omega, s * k, G),
                                 G);
            acc += val;
          }
          next[seg + k * m + q] = acc;
        }
      }
    }
    vec.swap(next);
  }
}

template <typename type>
static void applyDIF32Radix4ToBlock(std::vector<typename type::RX>& vec,
                                    const typename type::RX& root32,
                                    const typename type::RX& G)
{
  using RX = typename type::RX;

  // Stage 1: radix-4, n=32, m=8
  {
    long n = 32, m = 8;
    RX root = root32;
    RX omega4 = rxPowerMod<type>(root, m, G);
    std::vector<RX> next(32);
    for (long seg = 0; seg < 32; seg += n) {
      for (long q : range(m)) {
        for (long k : range(4)) {
          RX acc;
          clear(acc);
          for (long s : range(4)) {
            RX val = rxMulMod<type>(vec[seg + q + m * s],
                                    rxPowerMod<type>(root, q * k, G),
                                    G);
            val = rxMulMod<type>(val,
                                 rxPowerMod<type>(omega4, s * k, G),
                                 G);
            acc += val;
          }
          next[seg + k * m + q] = acc;
        }
      }
    }
    vec.swap(next);
  }

  // Stage 2: radix-4, n=8, m=2
  {
    long n = 8, m = 2;
    RX root = rxPowerMod<type>(root32, 32 / n, G);
    RX omega4 = rxPowerMod<type>(root, m, G);
    std::vector<RX> next(32);
    for (long seg = 0; seg < 32; seg += n) {
      for (long q : range(m)) {
        for (long k : range(4)) {
          RX acc;
          clear(acc);
          for (long s : range(4)) {
            RX val = rxMulMod<type>(vec[seg + q + m * s],
                                    rxPowerMod<type>(root, q * k, G),
                                    G);
            val = rxMulMod<type>(val,
                                 rxPowerMod<type>(omega4, s * k, G),
                                 G);
            acc += val;
          }
          next[seg + k * m + q] = acc;
        }
      }
    }
    vec.swap(next);
  }

  // Stage 3: radix-2, n=2, m=1
  {
    long n = 2, m = 1;
    RX root = rxPowerMod<type>(root32, 32 / n, G);
    RX omega = rxPowerMod<type>(root, m, G);
    std::vector<RX> next(32);
    for (long seg = 0; seg < 32; seg += n) {
      for (long q : range(m)) {
        for (long k : range(2)) {
          RX acc;
          clear(acc);
          for (long s : range(2)) {
            RX val = rxMulMod<type>(vec[seg + q + m * s],
                                    rxPowerMod<type>(root, q * k, G),
                                    G);
            val = rxMulMod<type>(val,
                                 rxPowerMod<type>(omega, s * k, G),
                                 G);
            acc += val;
          }
          next[seg + k * m + q] = acc;
        }
      }
    }
    vec.swap(next);
  }
}

template <typename type>
static NTL::Mat<typename type::RX>
buildRaderBlockMixMatrix(const NTL::Vec<long>& reps,
                         long cofactor,
                         const typename type::RX& root32,
                         const typename type::RX& G)
{
  using RX = typename type::RX;
  const long D = 96;
  NTL::Mat<RX> A = zeroMatrix<type>(D);

  std::vector<std::vector<RX>> bhat(3, std::vector<RX>(32));
  for (long t : range(D)) {
    long block = t % 3;
    long idx = t % 32;
    bhat[block][idx] = RX(reps[t] * cofactor, 1) % G;
  }

  for (long block : range(3))
    applyDIF32ToBlock<type>(bhat[block], root32, G);

  for (long out_block : range(3)) {
    for (long in_block : range(3)) {
      long b_block = mcMod(out_block - in_block, 3);
      for (long k : range(32)) {
        long in = in_block * 32 + k;
        long out = out_block * 32 + k;
        A[in][out] = bhat[b_block][k];
      }
    }
  }

  return A;
}

template <typename type>
static NTL::Mat<typename type::RX>
buildRaderBlockMixRadix4Matrix(const NTL::Vec<long>& reps,
                               long cofactor,
                               const typename type::RX& root32,
                               const typename type::RX& G)
{
  using RX = typename type::RX;
  const long D = 96;
  NTL::Mat<RX> A = zeroMatrix<type>(D);

  std::vector<std::vector<RX>> bhat(3, std::vector<RX>(32));
  for (long t : range(D)) {
    long block = t % 3;
    long idx = t % 32;
    bhat[block][idx] = RX(reps[t] * cofactor, 1) % G;
  }

  for (long block : range(3))
    applyDIF32Radix4ToBlock<type>(bhat[block], root32, G);

  for (long out_block : range(3)) {
    for (long in_block : range(3)) {
      long b_block = mcMod(out_block - in_block, 3);
      for (long k : range(32)) {
        long in = in_block * 32 + k;
        long out = out_block * 32 + k;
        A[in][out] = bhat[b_block][k];
      }
    }
  }

  return A;
}

template <typename type>
class ThinStep2Rader97Exec : public MatMulExecBase
{
  PA_INJECT(type)

  const EncryptedArray& base_ea;
  long dim;
  std::vector<std::unique_ptr<SparseMatMul1DExec>> stages;

public:
  ThinStep2Rader97Exec(const EncryptedArray& _ea,
                       std::shared_ptr<CubeSignature> sig,
                       const NTL::Vec<long>& reps,
                       long _dim,
                       long cofactor) :
      base_ea(_ea), dim(_dim)
  {
    std::cerr << "ThinStep2Rader97Exec: constructing for dim=" << dim << "\n";
    const EncryptedArrayDerived<type>& ea UNUSED = _ea.getDerived(type());
    RBak bak;
    bak.save();
    _ea.getAlMod().restoreContext();

    const RX& G = ea.getG();
    long D = sig->getDim(dim);
    std::cerr << "  D=" << D << ", deg(G)=" << NTL::deg(G) << "\n";
    assertEq(D, 96l, "ThinStep2Rader97Exec expects D=96");
    assertEq(reps.length(), D, "ThinStep2Rader97Exec expects 96 reps");

    R root32_R;
    std::cerr << "  Finding primitive 32nd root of unity in F_p...\n";
    FindPrimitiveRoot(root32_R, 32);
    std::cerr << "  Found root32.\n";
    RX root32 = rxConst<type>(root32_R);

    addStage<type>(stages, base_ea, dim, buildRaderInputCRTMatrix<type>(reps));
    std::cerr << "  Stage 1 (RaderInputCRT) added, " << stages.size() << " stages\n";

    for (long n = 32; n >= 2; n /= 2) {
      RX root = rxPowerMod<type>(root32, 32 / n, G);
      addStage<type>(stages,
                     base_ea,
                     dim,
                     buildDIF32StageMatrix<type>(n, root, G));
      std::cerr << "  DIF32 stage (n=" << n << ") added, " << stages.size() << " stages\n";
    }

    addStage<type>(stages,
                   base_ea,
                   dim,
                   buildRaderBlockMixMatrix<type>(reps, cofactor, root32, G));
    std::cerr << "  BlockMix stage added, " << stages.size() << " stages\n";

    for (long n = 2; n <= 32; n *= 2) {
      RX root = rxPowerMod<type>(root32, 32 / n, G);
      addStage<type>(stages,
                     base_ea,
                     dim,
                     buildInverseDIF32StageMatrix<type>(n, root, G));
      std::cerr << "  InvDIF32 stage (n=" << n << ") added, " << stages.size() << " stages\n";
    }

    addStage<type>(stages, base_ea, dim, buildCRTOutputMatrix<type>());
    std::cerr << "  CRTOutput stage added, " << stages.size() << " stages total\n";
    std::cerr << "  ThinStep2Rader97Exec construction complete.\n";

    // Compute partial products and report diagonal counts for grouping analysis
    {
      const RX& G_local = ea.getG();
      long D_local = 96;

      // Rebuild stage matrices for product computation
      std::vector<NTL::Mat<RX>> stage_mats;
      stage_mats.push_back(buildRaderInputCRTMatrix<type>(reps));
      for (long n = 32; n >= 2; n /= 2) {
        RX root_local = rxPowerMod<type>(root32, 32 / n, G_local);
        stage_mats.push_back(buildDIF32StageMatrix<type>(n, root_local, G_local));
      }
      stage_mats.push_back(buildRaderBlockMixMatrix<type>(reps, cofactor, root32, G_local));
      for (long n = 2; n <= 32; n *= 2) {
        RX root_local = rxPowerMod<type>(root32, 32 / n, G_local);
        stage_mats.push_back(buildInverseDIF32StageMatrix<type>(n, root_local, G_local));
      }
      stage_mats.push_back(buildCRTOutputMatrix<type>());

      // Compute cumulative products from left
      std::vector<long> cum_diag_counts;
      NTL::Mat<RX> product = stage_mats[0];
      for (long s = 1; s < (long)stage_mats.size(); s++) {
        NTL::Mat<RX> tmp;
        tmp.SetDims(D_local, D_local);
        for (long r = 0; r < D_local; r++)
          for (long c = 0; c < D_local; c++) {
            RX sum_val;
            NTL::clear(sum_val);
            for (long k = 0; k < D_local; k++) {
              if (IsZero(product[r][k]) || IsZero(stage_mats[s][k][c]))
                continue;
              sum_val += rxMulMod<type>(product[r][k], stage_mats[s][k][c], G_local);
            }
            tmp[r][c] = sum_val;
          }
        product = tmp;

        // Count non-zero diagonals
        long nz = 0;
        for (long d = 0; d < D_local; d++) {
          bool has_nz = false;
          for (long i = 0; i < D_local && !has_nz; i++)
            if (!IsZero(product[i][(i+d) % D_local]))
              has_nz = true;
          if (has_nz) nz++;
        }
        cum_diag_counts.push_back(nz);
        std::cerr << "  Product of stages 0.." << s << ": " << nz << " non-zero diagonals\n";
      }
    }
  }

  const EncryptedArray& getEA() const override { return base_ea; }

  void upgrade() override
  {
    for (auto& stage : stages)
      stage->upgrade();
  }

  void mul(Ctxt& ctxt) const override
  {
    HELIB_NTIMER_START(mul_ThinStep2Rader97Exec);
    for (const auto& stage : stages)
      stage->mul(ctxt);
  }
};

// ============================================================
// Mixed-Radix-4 DIF: 32 = 4 × 4 × 2 (9 stages instead of 13)
// ============================================================

template <typename type>
class ThinStep2Radix4Exec : public MatMulExecBase
{
  PA_INJECT(type)

  const EncryptedArray& base_ea;
  long dim;
  std::vector<std::unique_ptr<SparseMatMul1DExec>> stages;

public:
  ThinStep2Radix4Exec(const EncryptedArray& _ea,
                      std::shared_ptr<CubeSignature> sig,
                      const NTL::Vec<long>& reps,
                      long _dim,
                      long cofactor) :
      base_ea(_ea), dim(_dim)
  {
    std::cerr << "ThinStep2Radix4Exec: constructing for dim=" << dim << "\n";
    const EncryptedArrayDerived<type>& ea UNUSED = _ea.getDerived(type());
    RBak bak;
    bak.save();
    _ea.getAlMod().restoreContext();

    const RX& G = ea.getG();
    long D = sig->getDim(dim);
    assertEq(D, 96l, "ThinStep2Radix4Exec expects D=96");
    assertEq(reps.length(), D, "ThinStep2Radix4Exec expects 96 reps");

    R root32_R;
    FindPrimitiveRoot(root32_R, 32);
    RX root32 = rxConst<type>(root32_R);

    // Stage 1: RaderInputCRT
    addStage<type>(stages, base_ea, dim, buildRaderInputCRTMatrix<type>(reps));
    std::cerr << "  Stage 1 (RaderInputCRT) added\n";

    // Forward DIF: radix-4 (n=32), radix-4 (n=8), radix-2 (n=2)
    {
      RX root = root32;
      addStage<type>(stages, base_ea, dim,
                     buildDIF32Radix4StageMatrix<type>(32, root, G));
      std::cerr << "  Stage 2 (DIF radix-4, n=32) added\n";
    }
    {
      RX root = rxPowerMod<type>(root32, 32 / 8, G);
      addStage<type>(stages, base_ea, dim,
                     buildDIF32Radix4StageMatrix<type>(8, root, G));
      std::cerr << "  Stage 3 (DIF radix-4, n=8) added\n";
    }
    {
      RX root = rxPowerMod<type>(root32, 32 / 2, G);
      addStage<type>(stages, base_ea, dim,
                     buildDIF32StageMatrix<type>(2, root, G));
      std::cerr << "  Stage 4 (DIF radix-2, n=2) added\n";
    }

    // Stage 5: BlockMix (using radix-4 DIF ordering for b-hat coefficients)
    addStage<type>(stages, base_ea, dim,
                   buildRaderBlockMixRadix4Matrix<type>(reps, cofactor, root32, G));
    std::cerr << "  Stage 5 (BlockMix-Radix4) added\n";

    // Inverse DIF: radix-2 (n=2), radix-4 (n=8), radix-4 (n=32)
    {
      RX root = rxPowerMod<type>(root32, 32 / 2, G);
      addStage<type>(stages, base_ea, dim,
                     buildInverseDIF32StageMatrix<type>(2, root, G));
      std::cerr << "  Stage 6 (InvDIF radix-2, n=2) added\n";
    }
    {
      RX root = rxPowerMod<type>(root32, 32 / 8, G);
      addStage<type>(stages, base_ea, dim,
                     buildInverseDIF32Radix4StageMatrix<type>(8, root, G));
      std::cerr << "  Stage 7 (InvDIF radix-4, n=8) added\n";
    }
    {
      RX root = root32;
      addStage<type>(stages, base_ea, dim,
                     buildInverseDIF32Radix4StageMatrix<type>(32, root, G));
      std::cerr << "  Stage 8 (InvDIF radix-4, n=32) added\n";
    }

    // Stage 9: CRTOutput
    addStage<type>(stages, base_ea, dim, buildCRTOutputMatrix<type>());
    std::cerr << "  Stage 9 (CRTOutput) added, " << stages.size()
              << " stages total\n";
    std::cerr << "  ThinStep2Radix4Exec construction complete.\n";
  }

  const EncryptedArray& getEA() const override { return base_ea; }

  void upgrade() override
  {
    for (auto& stage : stages)
      stage->upgrade();
  }

  void mul(Ctxt& ctxt) const override
  {
    HELIB_NTIMER_START(mul_ThinStep2Radix4Exec);
    for (const auto& stage : stages)
      stage->mul(ctxt);
  }
};

// ============================================================
// Pure Power-of-2 Radix-4: D = 2^k (e.g., D=256=4^4)
// Rader converts (D+1)-point DFT to D-point cyclic convolution.
// DIF is a standard radix-4 FFT over the full D elements.
// ============================================================

template <typename type>
static NTL::Mat<typename type::RX>
buildDIFRadix4StageMatrix_Po2(long D, long n,
                              const typename type::RX& root,
                              const typename type::RX& G)
{
  using RX = typename type::RX;
  long m = n / 4;
  NTL::Mat<RX> A = zeroMatrix<type>(D);
  RX omega4 = rxPowerMod<type>(root, m, G);

  for (long seg = 0; seg < D; seg += n) {
    for (long q : range(m)) {
      for (long s : range(4)) {
        long in = seg + q + m * s;
        for (long k : range(4)) {
          long out = seg + k * m + q;
          RX val = rxMulMod<type>(rxPowerMod<type>(root, q * k, G),
                                   rxPowerMod<type>(omega4, s * k, G),
                                   G);
          A[in][out] += val;
        }
      }
    }
  }
  return A;
}

template <typename type>
static NTL::Mat<typename type::RX>
buildInvDIFRadix4StageMatrix_Po2(long D, long n,
                                 const typename type::RX& root,
                                 const typename type::RX& G)
{
  using R = typename type::R;
  using RX = typename type::RX;
  long m = n / 4;
  NTL::Mat<RX> A = zeroMatrix<type>(D);
  RX omega4 = rxPowerMod<type>(root, m, G);
  R inv4_R = inv(R(4));
  RX inv4 = rxConst<type>(inv4_R);

  for (long seg = 0; seg < D; seg += n) {
    for (long q : range(m)) {
      for (long s : range(4)) {
        long out = seg + q + m * s;
        for (long k : range(4)) {
          long in = seg + k * m + q;
          RX val = rxPowerMod<type>(root, mcMod(-q * k, n), G);
          val = rxMulMod<type>(val,
                               rxPowerMod<type>(omega4, mcMod(-s * k, 4), G),
                               G);
          val = rxMulMod<type>(val, inv4, G);
          A[in][out] += val;
        }
      }
    }
  }
  return A;
}

template <typename type>
static NTL::Mat<typename type::RX>
buildDIFRadix2StageMatrix_Po2(long D, long n,
                              const typename type::RX& root,
                              const typename type::RX& G)
{
  using RX = typename type::RX;
  long m = n / 2;
  NTL::Mat<RX> A = zeroMatrix<type>(D);
  RX omega = rxPowerMod<type>(root, m, G);

  for (long seg = 0; seg < D; seg += n) {
    for (long q : range(m)) {
      for (long s : range(2)) {
        long in = seg + q + m * s;
        for (long k : range(2)) {
          long out = seg + k * m + q;
          RX val = rxMulMod<type>(rxPowerMod<type>(root, q * k, G),
                                   rxPowerMod<type>(omega, s * k, G),
                                   G);
          A[in][out] += val;
        }
      }
    }
  }
  return A;
}

template <typename type>
static NTL::Mat<typename type::RX>
buildInvDIFRadix2StageMatrix_Po2(long D, long n,
                                 const typename type::RX& root,
                                 const typename type::RX& G)
{
  using R = typename type::R;
  using RX = typename type::RX;
  long m = n / 2;
  NTL::Mat<RX> A = zeroMatrix<type>(D);
  RX omega = rxPowerMod<type>(root, m, G);
  R inv2_R = inv(R(2));
  RX inv2 = rxConst<type>(inv2_R);

  for (long seg = 0; seg < D; seg += n) {
    for (long q : range(m)) {
      for (long s : range(2)) {
        long out = seg + q + m * s;
        for (long k : range(2)) {
          long in = seg + k * m + q;
          RX val = rxPowerMod<type>(root, mcMod(-q * k, n), G);
          val = rxMulMod<type>(val,
                               rxPowerMod<type>(omega, mcMod(-s * k, 2), G),
                               G);
          val = rxMulMod<type>(val, inv2, G);
          A[in][out] += val;
        }
      }
    }
  }
  return A;
}

template <typename type>
static NTL::Mat<typename type::RX>
buildRaderInputCRTMatrix_Po2(long D, const NTL::Vec<long>& reps)
{
  using RX = typename type::RX;
  NTL::Mat<RX> A = zeroMatrix<type>(D);
  RX one = rxConst<type>(1);
  RX minus_one = rxConst<type>(-1);

  for (long t : range(D)) {
    A[0][t] += minus_one;
    long coeff_idx = reps[mcMod(-t, D)];
    if (coeff_idx < D)
      A[coeff_idx][t] += one;
  }
  return A;
}

template <typename type>
static NTL::Mat<typename type::RX>
buildCRTOutputMatrix_Po2(long D)
{
  using RX = typename type::RX;
  NTL::Mat<RX> A = zeroMatrix<type>(D);
  RX one = rxConst<type>(1);
  for (long t : range(D))
    A[t][t] = one;
  return A;
}

template <typename type>
static NTL::Mat<typename type::RX>
buildRaderPointwiseMulMatrix_Po2(long D,
                                 const NTL::Vec<long>& reps,
                                 long cofactor,
                                 const typename type::RX& rootD,
                                 const typename type::RX& G)
{
  using RX = typename type::RX;
  NTL::Mat<RX> A = zeroMatrix<type>(D);

  // Compute b-hat = DIF(b) where b[k] = zeta^{reps[k]*cofactor}
  std::vector<RX> bvec(D);
  for (long t : range(D))
    bvec[t] = RX(reps[t] * cofactor, 1) % G;

  // Apply radix-4 DIF to bvec
  long logD = 0;
  { long tmp = D; while (tmp > 1) { tmp /= 2; logD++; } }
  long num_r4 = logD / 2;
  long has_r2 = logD % 2;

  for (long stage = 0; stage < num_r4; stage++) {
    long n = D >> (2 * stage);
    RX root = rxPowerMod<type>(rootD, D / n, G);
    RX omega4 = rxPowerMod<type>(root, n / 4, G);
    long m = n / 4;
    std::vector<RX> next(D);
    for (long seg = 0; seg < D; seg += n) {
      for (long q : range(m)) {
        for (long k : range(4)) {
          RX acc;
          clear(acc);
          for (long s : range(4)) {
            RX val = rxMulMod<type>(bvec[seg + q + m * s],
                                    rxPowerMod<type>(root, q * k, G), G);
            val = rxMulMod<type>(val, rxPowerMod<type>(omega4, s * k, G), G);
            acc += val;
          }
          next[seg + k * m + q] = acc;
        }
      }
    }
    bvec.swap(next);
  }
  if (has_r2) {
    long n = D >> (2 * num_r4);
    RX root = rxPowerMod<type>(rootD, D / n, G);
    RX omega = rxPowerMod<type>(root, n / 2, G);
    long m = n / 2;
    std::vector<RX> next(D);
    for (long seg = 0; seg < D; seg += n) {
      for (long q : range(m)) {
        for (long k : range(2)) {
          RX acc;
          clear(acc);
          for (long s : range(2)) {
            RX val = rxMulMod<type>(bvec[seg + q + m * s],
                                    rxPowerMod<type>(root, q * k, G), G);
            val = rxMulMod<type>(val, rxPowerMod<type>(omega, s * k, G), G);
            acc += val;
          }
          next[seg + k * m + q] = acc;
        }
      }
    }
    bvec.swap(next);
  }

  // Pointwise: diagonal matrix with bhat values
  for (long k : range(D))
    A[k][k] = bvec[k];

  return A;
}

template <typename type>
class ThinStep2Radix4Po2Exec : public MatMulExecBase
{
  PA_INJECT(type)

  const EncryptedArray& base_ea;
  long dim;
  // Hybrid: BSGS for dense stages, Sparse for DIF stages
  std::unique_ptr<MatMul1DExec> crt_in_bsgs;
  std::vector<std::unique_ptr<SparseMatMul1DExec>> dif_stages;
  std::unique_ptr<SparseMatMul1DExec> pointwise_stage;
  std::vector<std::unique_ptr<SparseMatMul1DExec>> inv_dif_stages;
  // CRT_out is identity for Po2, skip it

public:
  ThinStep2Radix4Po2Exec(const EncryptedArray& _ea,
                          std::shared_ptr<CubeSignature> sig,
                          const NTL::Vec<long>& reps,
                          long _dim,
                          long cofactor) :
      base_ea(_ea), dim(_dim)
  {
    const EncryptedArrayDerived<type>& ea UNUSED = _ea.getDerived(type());
    RBak bak;
    bak.save();
    _ea.getAlMod().restoreContext();

    const RX& G = ea.getG();
    long D = sig->getDim(dim);
    long logD = 0;
    { long tmp = D; while (tmp > 1) { tmp /= 2; logD++; } }
    assertTrue((1L << logD) == D, "ThinStep2Radix4Po2Exec: D must be power of 2");
    long num_r4 = logD / 2;
    long has_r2 = logD % 2;

    std::cerr << "ThinStep2Radix4Po2Exec(hybrid): D=" << D << ", logD=" << logD
              << ", radix-4 stages=" << num_r4
              << ", radix-2 stages=" << has_r2 << "\n";

    R rootD_R;
    FindPrimitiveRoot(rootD_R, D);
    RX rootD = rxConst<type>(rootD_R);

    // CRT_in: use BSGS (MatMul1DExec) for the dense D-diagonal matrix
    {
      NTL::Mat<RX> A = buildRaderInputCRTMatrix_Po2<type>(D, reps);
      Dense1DMatrix<type> mat(base_ea, dim, A);
      crt_in_bsgs.reset(new MatMul1DExec(mat, false));
      std::cerr << "  CRT_in (BSGS, " << D << " diags) added\n";
    }

    // Forward DIF: sparse stages (7 diags each)
    for (long stage = 0; stage < num_r4; stage++) {
      long n = D >> (2 * stage);
      RX root = rxPowerMod<type>(rootD, D / n, G);
      NTL::Mat<RX> A = buildDIFRadix4StageMatrix_Po2<type>(D, n, root, G);
      Dense1DMatrix<type> mat(base_ea, dim, A);
      dif_stages.emplace_back(new SparseMatMul1DExec(mat, diagonalOffsets<type>(A)));
      std::cerr << "  DIF radix-4 (n=" << n << ", "
                << diagonalOffsets<type>(A).size() << " diags) added\n";
    }
    if (has_r2) {
      long n = D >> (2 * num_r4);
      RX root = rxPowerMod<type>(rootD, D / n, G);
      NTL::Mat<RX> A = buildDIFRadix2StageMatrix_Po2<type>(D, n, root, G);
      Dense1DMatrix<type> mat(base_ea, dim, A);
      dif_stages.emplace_back(new SparseMatMul1DExec(mat, diagonalOffsets<type>(A)));
      std::cerr << "  DIF radix-2 (n=" << n << ", "
                << diagonalOffsets<type>(A).size() << " diags) added\n";
    }

    // Pointwise multiply: diagonal matrix (1 diag)
    {
      NTL::Mat<RX> A = buildRaderPointwiseMulMatrix_Po2<type>(D, reps, cofactor, rootD, G);
      Dense1DMatrix<type> mat(base_ea, dim, A);
      pointwise_stage.reset(new SparseMatMul1DExec(mat, diagonalOffsets<type>(A)));
      std::cerr << "  Pointwise (" << diagonalOffsets<type>(A).size() << " diags) added\n";
    }

    // Inverse DIF: sparse stages
    if (has_r2) {
      long n = D >> (2 * num_r4);
      RX root = rxPowerMod<type>(rootD, D / n, G);
      NTL::Mat<RX> A = buildInvDIFRadix2StageMatrix_Po2<type>(D, n, root, G);
      Dense1DMatrix<type> mat(base_ea, dim, A);
      inv_dif_stages.emplace_back(new SparseMatMul1DExec(mat, diagonalOffsets<type>(A)));
      std::cerr << "  InvDIF radix-2 (n=" << n << ") added\n";
    }
    for (long stage = num_r4 - 1; stage >= 0; stage--) {
      long n = D >> (2 * stage);
      RX root = rxPowerMod<type>(rootD, D / n, G);
      NTL::Mat<RX> A = buildInvDIFRadix4StageMatrix_Po2<type>(D, n, root, G);
      Dense1DMatrix<type> mat(base_ea, dim, A);
      inv_dif_stages.emplace_back(new SparseMatMul1DExec(mat, diagonalOffsets<type>(A)));
      std::cerr << "  InvDIF radix-4 (n=" << n << ") added\n";
    }

    std::cerr << "  Hybrid executor: 1 BSGS + "
              << dif_stages.size() + 1 + inv_dif_stages.size()
              << " sparse stages\n";
  }

  const EncryptedArray& getEA() const override { return base_ea; }

  void upgrade() override
  {
    crt_in_bsgs->upgrade();
    for (auto& s : dif_stages) s->upgrade();
    pointwise_stage->upgrade();
    for (auto& s : inv_dif_stages) s->upgrade();
  }

  void mul(Ctxt& ctxt) const override
  {
    HELIB_NTIMER_START(mul_ThinStep2Radix4Po2Exec);
    // CRT_in via BSGS (fast, sqrt(D) hoistings)
    crt_in_bsgs->mul(ctxt);
    // DIF stages via sparse (1 hoisting each, 7 diags)
    for (const auto& s : dif_stages) s->mul(ctxt);
    // Pointwise (1 hoisting, 1 diag)
    pointwise_stage->mul(ctxt);
    // InvDIF stages via sparse
    for (const auto& s : inv_dif_stages) s->mul(ctxt);
  }
};

static MatMulExecBase* buildThinStep2Radix4Po2Exec(const EncryptedArray& ea,
                                                   std::shared_ptr<CubeSignature> sig,
                                                   const NTL::Vec<long>& reps,
                                                   long dim,
                                                   long cofactor,
                                                   bool invert,
                                                   bool inflate)
{
  if (!thinStep2Radix4Enabled() || invert || inflate)
    return nullptr;
  if (ea.getTag() != PA_zz_p_tag)
    return nullptr;
  long D = sig->getDim(dim);
  if (D != (long)reps.length())
    return nullptr;
  // Check D is a power of 2 and D >= 16
  long logD = 0;
  { long tmp = D; while (tmp > 1) { tmp /= 2; logD++; } }
  if ((1L << logD) != D || D < 16)
    return nullptr;
  // Check that D-th root of unity exists in F_p (i.e., D | p-1)
  long p = ea.getAlMod().getZMStar().getP();
  if ((p - 1) % D != 0)
    return nullptr;

  return new ThinStep2Radix4Po2Exec<PA_zz_p>(ea, sig, reps, dim, cofactor);
}

// ============================================================
// Tensor Product Basis: I_6 ⊗ T_16 decomposition
// ============================================================

static bool tensorBasisEnabled()
{
  const char* flag = std::getenv("HELIB_TENSOR_BASIS");
  return flag != nullptr && std::atol(flag) != 0;
}

template <typename type>
class ThinStep2TensorExec : public MatMulExecBase
{
  PA_INJECT(type)

  const EncryptedArray& base_ea;
  long dim;
  std::unique_ptr<SparseMatMul1DExec> fallback_stage;

public:
  ThinStep2TensorExec(const EncryptedArray& _ea,
                      std::shared_ptr<CubeSignature> sig,
                      UNUSED const NTL::Vec<long>& reps,
                      long _dim,
                      UNUSED long cofactor,
                      bool invert) :
      base_ea(_ea), dim(_dim)
  {
    const EncryptedArrayDerived<type>& ea UNUSED = _ea.getDerived(type());
    RBak bak;
    bak.save();
    _ea.getAlMod().restoreContext();

    long D = sig->getDim(dim);
    long p = _ea.getAlMod().getZMStar().getP();

    std::cerr << "ThinStep2TensorExec(stride-BSGS): D=" << D << ", p=" << p
              << ", invert=" << invert << "\n";

    const auto& T = invert ? tensor_basis::T16_INV : tensor_basis::T16_FWD;

    // Build the 96x96 matrix with 16 non-zero diagonals at stride 6
    NTL::Mat<RX> A;
    A.SetDims(D, D);
    for (long d = 0; d < 16; d++) {
      for (long i = 0; i < D; i++) {
        long r = i / 6;
        long k = i % 6;
        long s = (r + d) % 16;
        long j = 6 * s + k;
        long val = ((long)T[r][s]) % p;
        if (val != 0)
          A[i][j] = rxConst<type>(val);
      }
    }

    // Use MatMul1DExec with the sparse matrix (it will use BSGS internally)
    // The key: even though D=96, only 16 diagonals are non-zero.
    // MatMul1DExec skips zero diagonals, so effective work = BSGS on 16 active diags.
    Dense1DMatrix<type> mat(base_ea, dim, A);
    fallback_stage.reset(new SparseMatMul1DExec(mat,
      std::vector<long>{0, 6, 12, 18, 24, 30, 36, 42, 48, 54, 60, 66, 72, 78, 84, 90}));

    std::cerr << "ThinStep2TensorExec: 16 diagonals at stride 6\n";
  }

  const EncryptedArray& getEA() const override { return base_ea; }

  void upgrade() override
  {
    if (fallback_stage) fallback_stage->upgrade();
  }

  void mul(Ctxt& ctxt) const override
  {
    HELIB_NTIMER_START(mul_ThinStep2TensorExec);
    if (fallback_stage) fallback_stage->mul(ctxt);
  }
};

static MatMulExecBase* buildThinStep2TensorExec(const EncryptedArray& ea,
                                                std::shared_ptr<CubeSignature> sig,
                                                const NTL::Vec<long>& reps,
                                                long dim,
                                                long cofactor,
                                                bool invert,
                                                bool inflate)
{
  if (!tensorBasisEnabled() || inflate)
    return nullptr;
  if (ea.getTag() != PA_zz_p_tag || sig->getDim(dim) != 96 ||
      reps.length() != 96)
    return nullptr;
  if (ea.getAlMod().getZMStar().getP() != 65537)
    return nullptr;  // Only for p=65537 (precomputed coefficients)

  return new ThinStep2TensorExec<PA_zz_p>(ea, sig, reps, dim, cofactor, invert);
}

// ============================================================

static MatMulExecBase* buildThinStep2Rader97Exec(const EncryptedArray& ea,
                                                 std::shared_ptr<CubeSignature> sig,
                                                 const NTL::Vec<long>& reps,
                                                 long dim,
                                                 long cofactor,
                                                 bool invert,
                                                 bool inflate)
{
  if (!thinStep2Rader97Enabled() || invert || inflate)
    return nullptr;
  if (ea.getTag() != PA_zz_p_tag || sig->getDim(dim) != 96 ||
      reps.length() != 96)
    return nullptr;

  return new ThinStep2Rader97Exec<PA_zz_p>(ea, sig, reps, dim, cofactor);
}

static MatMulExecBase* buildThinStep2Radix4Exec(const EncryptedArray& ea,
                                                std::shared_ptr<CubeSignature> sig,
                                                const NTL::Vec<long>& reps,
                                                long dim,
                                                long cofactor,
                                                bool invert,
                                                bool inflate)
{
  if (!thinStep2Radix4Enabled() || invert || inflate)
    return nullptr;
  if (ea.getTag() != PA_zz_p_tag || sig->getDim(dim) != 96 ||
      reps.length() != 96)
    return nullptr;

  return new ThinStep2Radix4Exec<PA_zz_p>(ea, sig, reps, dim, cofactor);
}

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
    const EncryptedArrayDerived<type>& ea UNUSED = _ea.getDerived(type());
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
