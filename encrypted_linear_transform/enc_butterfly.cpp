/* enc_butterfly.cpp
 *
 * Encrypted prototype of the O(log D) mixed-radix/Rader butterfly evaluation
 * of the CoeffToSlot (ThinEvalMap Step2) dim-0 block for HElib parameter
 * Set E: p=65537, m=50731=97*523, D=96, slot ring F_p[X]/G, deg G = 18.
 *
 * Pipeline (verified in plaintext by derive_butterfly.py):
 *   w = INTT3(INTT2(INTT1( Dmid . NTT3(NTT2(NTT1(v))) )))     [6 stages]
 *   out = w with slot (dim0==48) replaced by DC = sum_j v[j]  [Rader DC fixup]
 *   True Step2 output: y[i] = out[dlog_5(i)] (i>=1), y[0]=out[48]
 *   i.e. encrypted output equals the true block output up to the fixed,
 *   known Rader permutation (absorbable into downstream linear maps).
 *
 * Also runs the monolithic BSGS MatMul1D of the same 96x96 block as baseline.
 *
 * Usage: enc_butterfly [bits=600] [allks=1] [baseline=1]
 */
#include <helib/helib.h>
#include <helib/matmul.h>
#include <helib/debugging.h>
#include <NTL/ZZX.h>
#include <NTL/lzz_pE.h>
#include <fstream>
#include <sstream>
#include <iostream>
#include <chrono>
#include <vector>

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
  // diags[offidx][slot_i] = vector<long> of d coeffs
  vector<vector<vector<long>>> diags;
};

int main(int argc, char** argv)
{
  long bits = 600, allks = 1, baseline = 1;
  if (argc > 1) bits = atol(argv[1]);
  if (argc > 2) allks = atol(argv[2]);
  if (argc > 3) baseline = atol(argv[3]);

  const long m = 50731, p = 65537, r = 1, c = 3;
  const long D = 96;
  const string proto_dir =
      "/home/user/experiments/order4bgv/tower_linear_transform/"
      "encrypted_prototype/";
  const string dump_path =
      "/home/user/order4-artifact/"
      "step2_matrix_dim0_sz96.txt";

  cout << "=== encrypted butterfly prototype: m=" << m << " p=" << p
       << " bits=" << bits << " allks=" << allks << " ===" << endl;

  // ---------- load stage data ----------
  ifstream sf(proto_dir + "stages.txt");
  if (!sf) { cerr << "cannot open stages.txt" << endl; return 1; }
  long fp, fD, fd;
  sf >> fp >> fD >> fd;
  if (fp != p || fD != D) { cerr << "stages.txt mismatch" << endl; return 1; }
  const long dd = fd; // 18
  vector<long> Gcoeffs(dd + 1);
  for (auto& x : Gcoeffs) sf >> x;
  long nstages;
  sf >> nstages;
  vector<StageData> stages(nstages);
  long total_rot = 0, total_mul = 0;
  for (auto& st : stages) {
    long noff; sf >> noff;
    st.offsets.resize(noff);
    st.diags.resize(noff);
    for (long oi = 0; oi < noff; oi++) {
      sf >> st.offsets[oi];
      st.diags[oi].resize(D);
      for (long i = 0; i < D; i++) {
        st.diags[oi][i].resize(dd);
        for (long k = 0; k < dd; k++) sf >> st.diags[oi][i][k];
      }
      if (st.offsets[oi] != 0) total_rot++;
      total_mul++;
    }
  }
  cout << "loaded " << nstages << " stages, nontrivial rotations in chain: "
       << total_rot << ", const-mults: " << total_mul << endl;

  // ---------- load dump matrix ----------
  vector<vector<vector<long>>> Smat(D, vector<vector<long>>(D));
  {
    ifstream df(dump_path);
    if (!df) { cerr << "cannot open dump" << endl; return 1; }
    string line;
    getline(df, line); // "96"
    for (long i = 0; i < D; i++) {
      getline(df, line);
      stringstream ss(line);
      string entry;
      for (long j = 0; j < D; j++) {
        getline(ss, entry, '|');
        vector<long> cs;
        for (size_t q = 0; q < entry.size(); q++)
          if (entry[q] == '[' || entry[q] == ']') entry[q] = ' ';
        stringstream es(entry);
        long v;
        while (es >> v) cs.push_back(v);
        cs.resize(dd, 0);
        Smat[i][j] = cs;
      }
    }
  }
  cout << "loaded dump matrix" << endl;

  // ---------- context ----------
  double t0 = now_sec();
  Context context = ContextBuilder<BGV>()
                        .m(m).p(p).r(r).bits(bits).c(c)
                        .gens(vector<long>{48117, 5239})
                        .ords(vector<long>{96, 29})
                        .build();
  cout << "context built in " << now_sec() - t0 << " s" << endl;
  cout << "securityLevel = " << context.securityLevel() << endl;
  const EncryptedArray& ea = context.getEA();
  long nslots = ea.size();
  long d_slot = context.getOrdP();
  cout << "nslots=" << nslots << " d=" << d_slot
       << " dim0 size=" << ea.sizeOfDimension(0)
       << " native=" << ea.nativeDimension(0)
       << " dim1 size=" << ea.sizeOfDimension(1)
       << " native=" << ea.nativeDimension(1) << endl;
  if (d_slot != dd) { cerr << "slot degree mismatch" << endl; return 1; }

  // verify HElib G == our G
  {
    const EncryptedArrayDerived<PA_zz_p>& ead = ea.getDerived(PA_zz_p());
    zz_pBak bak; bak.save();
    context.getAlMod().restoreContext();
    const zz_pX& G = ead.getG();
    bool okG = (deg(G) == dd);
    for (long k = 0; k <= dd && okG; k++)
      if (rep(coeff(G, k)) != Gcoeffs[k] % p) okG = false;
    cout << "HElib slot polynomial G matches stage data G: "
         << (okG ? "YES" : "NO") << endl;
    if (!okG) {
      cout << "  HElib G: ";
      for (long k = 0; k <= deg(G); k++) cout << rep(coeff(G, k)) << " ";
      cout << endl << "  file  G: ";
      for (auto& x : Gcoeffs) cout << x << " ";
      cout << endl;
      return 1;
    }
  }

  // ---------- keys ----------
  t0 = now_sec();
  SecKey secretKey(context);
  secretKey.GenSecKey();
  if (allks)
    add1DMatrices(secretKey);   // KS matrices for ALL 1D rotation amounts
  else
    addSome1DMatrices(secretKey);
  const PubKey& publicKey = secretKey;
  cout << "keygen+KS in " << now_sec() - t0 << " s" << endl;

  // ---------- plaintext input ----------
  SetSeed(ZZ(20260722));
  vector<ZZX> vslots(nslots);
  const PAlgebra& zms = context.getZMStar();
  for (long k = 0; k < nslots; k++) {
    ZZX f;
    for (long q = 0; q < dd; q++) SetCoeff(f, q, RandomBnd(p));
    vslots[k] = f;
  }

  // slot index <-> (i0,i1) coordinates
  vector<long> coord0(nslots), coord1(nslots);
  for (long k = 0; k < nslots; k++) {
    coord0[k] = zms.coordinate(0, k);
    coord1[k] = zms.coordinate(1, k);
  }

  // ---------- encrypt ----------
  Ctxt ctxt(publicKey);
  ea.encrypt(ctxt, publicKey, vslots);
  cout << "input capacity: " << ctxt.capacity() << " bits" << endl;

  // ---------- rotation direction probe ----------
  // encrypt slot value = coord0, rotate1D(.,0,1), see where content lands
  long sgn = 0;
  {
    vector<ZZX> probe(nslots);
    for (long k = 0; k < nslots; k++) probe[k] = ZZX(coord0[k]);
    Ctxt pc(publicKey);
    ea.encrypt(pc, publicKey, probe);
    ea.rotate1D(pc, 0, 1);
    vector<ZZX> dec;
    ea.decrypt(pc, secretKey, dec);
    // find a slot with coord0==5, check its value
    for (long k = 0; k < nslots; k++)
      if (coord0[k] == 5) {
        long val = to_long(coeff(dec[k], 0));
        if (val == 4) sgn = -1;      // slot i now holds old i-1 => rot(+k): i <- i-k
        else if (val == 6) sgn = +1; // slot i holds old i+1
        break;
      }
    if (sgn == 0) { cerr << "direction probe failed" << endl; return 1; }
    // to obtain v[(i+o)%D] at slot i we need rotate1D by amt = sgn_amt:
    // if rot(+k) gives slot i = old[i-k], then rot(-o) gives old[i+o] -> amt=-o
    cout << "rotate1D(+1) => slot i holds old[i" << (sgn == -1 ? "-" : "+")
         << "1]; using amt = " << (sgn == -1 ? "-offset" : "+offset") << endl;
  }
  long rotamt_factor = (sgn == -1) ? -1 : +1; // amt = rotamt_factor * (-o)... set below
  // we need old[(i+o)]: if rot(+k): i<-old[i-k], choose k=-o. So amt = -o.
  // if rot(+k): i<-old[i+k], choose k=+o. So amt = +o.
  long amt_mult = (sgn == -1) ? -1 : +1;

  // ---------- pre-encode all stage constants ----------
  t0 = now_sec();
  vector<vector<ZZX>> encoded(nstages);
  for (long s = 0; s < nstages; s++) {
    encoded[s].resize(stages[s].offsets.size());
    for (size_t oi = 0; oi < stages[s].offsets.size(); oi++) {
      vector<ZZX> cvec(nslots);
      for (long k = 0; k < nslots; k++) {
        ZZX f;
        for (long q = 0; q < dd; q++)
          SetCoeff(f, q, stages[s].diags[oi][coord0[k]][q]);
        cvec[k] = f;
      }
      ea.encode(encoded[s][oi], cvec);
    }
  }
  // masks for DC insertion at dim0 slot 48
  ZZX mask48, masknot48;
  {
    vector<long> m1(nslots), m2(nslots);
    for (long k = 0; k < nslots; k++) {
      m1[k] = (coord0[k] == 48) ? 1 : 0;
      m2[k] = 1 - m1[k];
    }
    ea.encode(mask48, m1);
    ea.encode(masknot48, m2);
  }
  cout << "constants encoded in " << now_sec() - t0 << " s" << endl;

  // ---------- encrypted butterfly chain ----------
  setTimersOn();
  long rot_count = 0;
  double t_chain0 = now_sec();
  Ctxt x = ctxt; // working copy
  Ctxt dc_in = ctxt; // for DC path (uses original input)
  for (long s = 0; s < nstages; s++) {
    double ts = now_sec();
    Ctxt acc(ZeroCtxtLike, x);
    for (size_t oi = 0; oi < stages[s].offsets.size(); oi++) {
      long o = stages[s].offsets[oi];
      Ctxt tmp = x;
      if (o != 0) { ea.rotate1D(tmp, 0, amt_mult * o); rot_count++; }
      tmp.multByConstant(encoded[s][oi]);
      acc += tmp;
    }
    x = acc;
    cout << "  stage " << s + 1 << ": " << now_sec() - ts << " s, capacity "
         << x.capacity() << " bits" << endl;
  }
  // DC path: total sum along dim 0 (7 rotations, amounts 32,64,16,8,4,2,1)
  {
    double ts = now_sec();
    Ctxt A = dc_in;
    Ctxt tmp = dc_in;
    ea.rotate1D(tmp, 0, 32); rot_count++;
    A += tmp;
    tmp = dc_in;
    ea.rotate1D(tmp, 0, 64); rot_count++;
    A += tmp;
    for (long amt : {16, 8, 4, 2, 1}) {
      tmp = A;
      ea.rotate1D(tmp, 0, amt); rot_count++;
      A += tmp;
    }
    // insert DC at dim0 position 48
    x.multByConstant(masknot48);
    A.multByConstant(mask48);
    x += A;
    cout << "  DC fixup: " << now_sec() - ts << " s, capacity "
         << x.capacity() << " bits" << endl;
  }
  double t_chain = now_sec() - t_chain0;
  cout << "ENCRYPTED CHAIN: " << t_chain << " s total, " << rot_count
       << " rotations, final capacity " << x.capacity() << " bits (input "
       << ctxt.capacity() << ")" << endl;

  // ---------- decrypt and verify ----------
  vector<ZZX> out;
  ea.decrypt(x, secretKey, out);

  // plaintext reference with NTL zz_pE
  zz_pBak bak; bak.save();
  zz_p::init(p);
  zz_pX Gx;
  for (long k = 0; k <= dd; k++) SetCoeff(Gx, k, Gcoeffs[k]);
  zz_pE::init(Gx);
  auto toE = [&](const vector<long>& cs) {
    zz_pX f;
    for (long q = 0; q < (long)cs.size(); q++) SetCoeff(f, q, cs[q]);
    return conv<zz_pE>(f);
  };
  auto zzxToE = [&](const ZZX& f) {
    zz_pX g;
    for (long q = 0; q <= deg(f); q++)
      SetCoeff(g, q, conv<zz_p>(coeff(f, q)));
    return conv<zz_pE>(g);
  };
  // dlog table: i = 5^t mod 97
  vector<long> pow5(96);
  {
    long gp = 1;
    for (long t = 0; t < 96; t++) { pow5[t] = gp; gp = (gp * 5) % 97; }
  }
  // for each dim1 slice: ref[i] = sum_j S[i][j] * v[(j,i1)]
  long n1 = ea.sizeOfDimension(1);
  vector<vector<long>> slot_at(D, vector<long>(n1, -1));
  for (long k = 0; k < nslots; k++) slot_at[coord0[k]][coord1[k]] = k;
  long mismatches = 0, checked = 0;
  for (long i1 = 0; i1 < n1; i1++) {
    vector<zz_pE> vin(D), ref(D);
    for (long j = 0; j < D; j++) vin[j] = zzxToE(vslots[slot_at[j][i1]]);
    for (long i = 0; i < D; i++) {
      zz_pE acc;
      clear(acc);
      for (long j = 0; j < D; j++) acc += toE(Smat[i][j]) * vin[j];
      ref[i] = acc;
    }
    // expected encrypted output: out[t] = ref[pow5[t]] for t != 48 (pow5[t]<96)
    //                            out[48] = ref[0]
    for (long t = 0; t < D; t++) {
      zz_pE expect = (t == 48) ? ref[0] : ref[pow5[t]];
      zz_pE got = zzxToE(out[slot_at[t][i1]]);
      checked++;
      if (expect != got) mismatches++;
    }
  }
  cout << "VERIFY vs dump matvec (up to fixed Rader permutation): "
       << (mismatches == 0 ? "PASS" : "FAIL") << "  (" << checked
       << " slots checked, " << mismatches << " mismatches)" << endl;

  // ---------- baseline: monolithic BSGS MatMul1D ----------
  if (baseline) {
    struct DumpMat : public MatMul1D_derived<PA_zz_p> {
      const EncryptedArray& ea_;
      const vector<vector<vector<long>>>& M;
      long dd_;
      bool transposed;
      DumpMat(const EncryptedArray& ea2,
              const vector<vector<vector<long>>>& M2,
              long dd2, bool tr)
          : ea_(ea2), M(M2), dd_(dd2), transposed(tr) {}
      const EncryptedArray& getEA() const override { return ea_; }
      long getDim() const override { return 0; }
      bool multipleTransforms() const override { return false; }
      bool get(zz_pX& out, long i, long j, long k) const override {
        const vector<long>& cs = transposed ? M[j][i] : M[i][j];
        clear(out);
        bool zero = true;
        for (long q = 0; q < dd_; q++)
          if (cs[q] % 65537) { SetCoeff(out, q, cs[q]); zero = false; }
        return zero;
      }
    };
    for (int tr = 0; tr <= 1; tr++) {
      DumpMat mat(ea, Smat, dd, tr != 0);
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
      // verify vs ref (recompute per slice, orientation as in our convention)
      long mm = 0;
      for (long i1 = 0; i1 < n1; i1++) {
        vector<zz_pE> vin(D), ref(D);
        for (long j = 0; j < D; j++) vin[j] = zzxToE(vslots[slot_at[j][i1]]);
        for (long i = 0; i < D; i++) {
          zz_pE acc;
          clear(acc);
          for (long j = 0; j < D; j++) acc += toE(Smat[i][j]) * vin[j];
          ref[i] = acc;
        }
        for (long i = 0; i < D; i++)
          if (ref[i] != zzxToE(bout[slot_at[i][i1]])) mm++;
      }
      cout << "BASELINE BSGS MatMul1D (transposed=" << tr << "): precomp "
           << tb_pre << " s, mul " << tb << " s, capacity after "
           << y.capacity() << " bits (input " << ctxt.capacity()
           << "), matvec match: " << (mm == 0 ? "YES" : "no") << " (" << mm
           << " mismatches)" << endl;
    }
  }

  printNamedTimer(cout, "AAA_rotate1D");
  cout << "done" << endl;
  return 0;
}
