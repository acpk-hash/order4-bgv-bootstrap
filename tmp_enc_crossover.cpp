/* enc_crossover.cpp
 *
 * Crossover microbenchmark: encrypted mixed-radix butterfly evaluation of the
 * length-D DFT over the slot field vs monolithic BSGS MatMul1D of the same
 * dense D x D matrix, on a ring whose hypercube has a dimension of size D.
 *
 * D is smooth (D = l-1, p == 1 mod lcm(l, l-1)) so zeta_D lives in F_p and the
 * butterfly is a PURE mixed-radix DIF chain (no Rader wrap, scalar constants).
 * Stage data + output permutation come from derive_crossover.py (plaintext-
 * verified: stage product == row-permuted DFT).
 *
 * Usage: enc_crossover <stagefile> <m> [bits=600] [baseline=1] [butterfly=1]
 */
#include <helib/helib.h>
#include <helib/matmul.h>
#include <helib/timing.h>
#include <NTL/ZZX.h>
#include <NTL/lzz_pX.h>
#include <fstream>
#include <sstream>
#include <iostream>
#include <chrono>
#include <vector>
#include <set>

using namespace helib;
using namespace NTL;
using namespace std;

static double now_sec()
{
  return chrono::duration<double>(
             chrono::steady_clock::now().time_since_epoch())
      .count();
}

struct StageData {
  vector<long> offsets;
  vector<vector<long>> diags; // diags[offidx][row i] scalar
};

int main(int argc, char** argv)
{
  if (argc < 3) {
    cerr << "usage: enc_crossover <stagefile> <m> [bits] [baseline] [butterfly]"
         << endl;
    return 1;
  }
  const string stagefile = argv[1];
  const long m = atol(argv[2]);
  long bits = (argc > 3) ? atol(argv[3]) : 600;
  long baseline = (argc > 4) ? atol(argv[4]) : 1;
  long butterfly = (argc > 5) ? atol(argv[5]) : 1;

  // ---------- load stage data ----------
  long p, D, zeta;
  long nstages;
  vector<StageData> stages;
  vector<long> perm;
  {
    ifstream sf(stagefile);
    if (!sf) { cerr << "cannot open " << stagefile << endl; return 1; }
    sf >> p >> D >> zeta >> nstages;
    stages.resize(nstages);
    for (auto& st : stages) {
      long noff; sf >> noff;
      st.offsets.resize(noff);
      st.diags.resize(noff);
      for (long oi = 0; oi < noff; oi++) {
        sf >> st.offsets[oi];
        st.diags[oi].resize(D);
        for (long i = 0; i < D; i++) sf >> st.diags[oi][i];
      }
    }
    perm.resize(D);
    for (long i = 0; i < D; i++) sf >> perm[i];
    if (!sf) { cerr << "stage file parse error" << endl; return 1; }
  }
  long chain_rot = 0, chain_mul = 0;
  for (auto& st : stages)
    for (long o : st.offsets) { if (o) chain_rot++; chain_mul++; }
  cout << "=== crossover benchmark: m=" << m << " p=" << p << " D=" << D
       << " bits=" << bits << " ===" << endl;
  cout << "butterfly chain: " << nstages << " stages, " << chain_rot
       << " nontrivial rotations, " << chain_mul << " const-mults" << endl;

  // zeta powers
  vector<long> zp(D);
  zp[0] = 1;
  for (long i = 1; i < D; i++) zp[i] = MulMod(zp[i - 1], zeta, p);

  // ---------- context ----------
  double t0 = now_sec();
  Context context =
      ContextBuilder<BGV>().m(m).p(p).r(1).bits(bits).c(3).build();
  cout << "context built in " << now_sec() - t0
       << " s, securityLevel=" << context.securityLevel() << endl;
  const EncryptedArray& ea = context.getEA();
  const PAlgebra& zms = context.getZMStar();
  long nslots = ea.size();
  long d_slot = context.getOrdP();
  long ndims = ea.dimension();
  long dD = -1;
  for (long i = 0; i < ndims; i++) {
    cout << "  dim " << i << ": size " << ea.sizeOfDimension(i)
         << " native=" << ea.nativeDimension(i) << endl;
    if (ea.sizeOfDimension(i) == D) dD = i;
  }
  cout << "nslots=" << nslots << " d=" << d_slot << " target dim=" << dD
       << endl;
  if (dD < 0) { cerr << "no dimension of size D" << endl; return 1; }
  if (!ea.nativeDimension(dD)) { cerr << "target dim not native" << endl; return 1; }

  // ---------- keys ----------
  t0 = now_sec();
  SecKey secretKey(context);
  secretKey.GenSecKey();
  addSome1DMatrices(secretKey);
  double t_ks_some = now_sec() - t0;
  // targeted KS matrices for all butterfly rotation amounts (both signs)
  t0 = now_sec();
  std::set<long> amts;
  for (auto& st : stages)
    for (long o : st.offsets)
      if (o) { amts.insert(o); amts.insert(D - o); }
  long gen_extra = 0;
  for (long a : amts) {
    long el = zms.genToPow(dD, a);
    if (!secretKey.haveKeySWmatrix(1, el)) {
      secretKey.GenKeySWmatrix(1, el);
      gen_extra++;
    }
  }
  secretKey.setKeySwitchMap();
  double t_ks_extra = now_sec() - t0;
  const PubKey& publicKey = secretKey;
  cout << "keygen: addSome1DMatrices " << t_ks_some << " s, +" << gen_extra
       << " targeted butterfly KS matrices " << t_ks_extra << " s" << endl;

  // ---------- plaintext input ----------
  SetSeed(ZZ(20260722));
  vector<ZZX> vslots(nslots);
  for (long k = 0; k < nslots; k++) {
    ZZX f;
    for (long q = 0; q < d_slot; q++) SetCoeff(f, q, RandomBnd(p));
    vslots[k] = f;
  }
  // coordinates: cD = coord in target dim, rest = flattened other coords
  vector<long> cD(nslots), rest(nslots);
  long nrest = nslots / D;
  for (long k = 0; k < nslots; k++) {
    cD[k] = zms.coordinate(dD, k);
    long r2 = 0;
    for (long i = 0; i < ndims; i++)
      if (i != dD) r2 = r2 * ea.sizeOfDimension(i) + zms.coordinate(i, k);
    rest[k] = r2;
  }
  vector<vector<long>> slot_at(D, vector<long>(nrest, -1));
  for (long k = 0; k < nslots; k++) slot_at[cD[k]][rest[k]] = k;
  for (long i = 0; i < D; i++)
    for (long r2 = 0; r2 < nrest; r2++)
      if (slot_at[i][r2] < 0) { cerr << "slot indexing hole" << endl; return 1; }

  // ---------- reference: ref[i] = sum_j W[i][j] * v[j], W[i][j]=zeta^(ij) ----
  // (coefficient-wise on the slot polynomials since W entries are scalars)
  double t_ref0 = now_sec();
  // refc[r2][i][q]
  vector<vector<vector<long>>> refc(
      nrest, vector<vector<long>>(D, vector<long>(d_slot, 0)));
  {
    vector<vector<long>> vc(D, vector<long>(d_slot));
    for (long r2 = 0; r2 < nrest; r2++) {
      for (long j = 0; j < D; j++) {
        const ZZX& f = vslots[slot_at[j][r2]];
        for (long q = 0; q < d_slot; q++)
          vc[j][q] = to_long(coeff(f, q) % p);
      }
      for (long i = 0; i < D; i++) {
        vector<long>& out = refc[r2][i];
        for (long j = 0; j < D; j++) {
          long w = zp[(i * j) % D];
          const vector<long>& vj = vc[j];
          for (long q = 0; q < d_slot; q++)
            out[q] = AddMod(out[q], MulMod(w, vj[q], p), p);
        }
      }
    }
  }
  cout << "plaintext reference computed in " << now_sec() - t_ref0 << " s"
       << endl;
  auto check_out = [&](const vector<ZZX>& out, bool permuted) {
    long mm = 0;
    for (long r2 = 0; r2 < nrest; r2++)
      for (long t = 0; t < D; t++) {
        long i = permuted ? perm[t] : t;
        const ZZX& f = out[slot_at[t][r2]];
        for (long q = 0; q < d_slot; q++) {
          long got = to_long(coeff(f, q) % p);
          long want = refc[r2][i][q];
          if (got != want) { mm++; break; }
        }
      }
    return mm;
  };

  // ---------- encrypt ----------
  Ctxt ctxt(publicKey);
  ea.encrypt(ctxt, publicKey, vslots);
  cout << "input capacity: " << ctxt.capacity() << " bits" << endl;

  // ---------- rotation direction probe ----------
  long amt_mult = 0;
  {
    vector<ZZX> probe(nslots);
    for (long k = 0; k < nslots; k++) probe[k] = ZZX(cD[k]);
    Ctxt pc(publicKey);
    ea.encrypt(pc, publicKey, probe);
    ea.rotate1D(pc, dD, 1);
    vector<ZZX> dec;
    ea.decrypt(pc, secretKey, dec);
    for (long k = 0; k < nslots; k++)
      if (cD[k] == 5) {
        long val = to_long(coeff(dec[k], 0));
        if (val == 4) amt_mult = -1;      // rot(+1): slot i <- old[i-1]
        else if (val == 6) amt_mult = +1; // rot(+1): slot i <- old[i+1]
        break;
      }
    if (amt_mult == 0) { cerr << "direction probe failed" << endl; return 1; }
    cout << "rotate1D(+1) => slot i holds old[i" << (amt_mult < 0 ? "-" : "+")
         << "1]; to fetch v[i+o] use amt = " << (amt_mult < 0 ? "-o" : "+o")
         << endl;
  }

  setTimersOn();

  // ---------- encrypted butterfly ----------
  if (butterfly) {
    // pre-encode constants
    t0 = now_sec();
    vector<vector<ZZX>> encoded(nstages);
    for (long s = 0; s < nstages; s++) {
      encoded[s].resize(stages[s].offsets.size());
      for (size_t oi = 0; oi < stages[s].offsets.size(); oi++) {
        vector<long> cvec(nslots);
        for (long k = 0; k < nslots; k++)
          cvec[k] = stages[s].diags[oi][cD[k]];
        ea.encode(encoded[s][oi], cvec);
      }
    }
    cout << "butterfly constants encoded in " << now_sec() - t0 << " s"
         << endl;
    resetAllTimers();
    long rot_count = 0;
    double t_chain0 = now_sec();
    Ctxt x = ctxt;
    for (long s = 0; s < nstages; s++) {
      double ts = now_sec();
      Ctxt acc(ZeroCtxtLike, x);
      for (size_t oi = 0; oi < stages[s].offsets.size(); oi++) {
        long o = stages[s].offsets[oi];
        Ctxt tmp = x;
        if (o != 0) {
          long amt = ((amt_mult * o) % D + D) % D;
          ea.rotate1D(tmp, dD, amt);
          rot_count++;
        }
        tmp.multByConstant(encoded[s][oi]);
        acc += tmp;
      }
      x = acc;
      cout << "  stage " << s + 1 << ": " << now_sec() - ts << " s, capacity "
           << x.capacity() << " bits" << endl;
    }
    double t_chain = now_sec() - t_chain0;
    cout << "BUTTERFLY: " << t_chain << " s, " << rot_count
         << " rotations, capacity " << ctxt.capacity() << " -> "
         << x.capacity() << " (consumed "
         << ctxt.capacity() - x.capacity() << " bits)" << endl;
    vector<ZZX> out;
    ea.decrypt(x, secretKey, out);
    long mm = check_out(out, true);
    cout << "BUTTERFLY VERIFY vs DFT matvec (up to fixed digit-reversal perm): "
         << (mm == 0 ? "PASS" : "FAIL") << " (" << mm
         << " mismatched slots of " << nslots << ")" << endl;
    cout << "--- butterfly timers ---" << endl;
    printNamedTimer(cout, "smartAutomorph");
    printNamedTimer(cout, "automorph");
    printNamedTimer(cout, "AAA_rotate1D");
  }

  // ---------- baseline: monolithic BSGS MatMul1D ----------
  if (baseline) {
    struct DFTMat : public MatMul1D_derived<PA_zz_p> {
      const EncryptedArray& ea_;
      const vector<long>& zp_;
      long p_, D_, dim_;
      bool transposed;
      DFTMat(const EncryptedArray& ea2, const vector<long>& zp2, long pp,
             long DD, long dim2, bool tr)
          : ea_(ea2), zp_(zp2), p_(pp), D_(DD), dim_(dim2), transposed(tr) {}
      const EncryptedArray& getEA() const override { return ea_; }
      long getDim() const override { return dim_; }
      bool multipleTransforms() const override { return false; }
      bool get(zz_pX& out, long i, long j, long k) const override {
        long a = transposed ? j : i, b = transposed ? i : j;
        long w = zp_[(a * b) % D_];
        clear(out);
        if (w == 0) return true;
        SetCoeff(out, 0, w);
        return false;
      }
    };
    for (int tr = 0; tr <= 1; tr++) {
      resetAllTimers();
      DFTMat mat(ea, zp, p, D, dD, tr != 0);
      double tb0 = now_sec();
      MatMul1DExec exec(mat, false);
      exec.upgrade();
      double tb_pre = now_sec() - tb0;
      Ctxt y = ctxt;
      tb0 = now_sec();
      exec.mul(y);
      double tb = now_sec() - tb0;
      vector<ZZX> bout;
      ea.decrypt(y, secretKey, bout);
      long mm = check_out(bout, false);
      cout << "BSGS MatMul1D (transposed=" << tr << "): precomp " << tb_pre
           << " s, mul " << tb << " s, capacity " << ctxt.capacity() << " -> "
           << y.capacity() << " (consumed " << ctxt.capacity() - y.capacity()
           << " bits), verify: " << (mm == 0 ? "PASS" : "no") << " (" << mm
           << " mismatched slots)" << endl;
      cout << "--- BSGS timers (transposed=" << tr << ") ---" << endl;
      printNamedTimer(cout, "smartAutomorph");
      printNamedTimer(cout, "automorph");
      printNamedTimer(cout, "AAA_rotate1D");
    }
  }
  cout << "done" << endl;
  return 0;
}
