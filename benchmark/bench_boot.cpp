// bench_boot.cpp
//
// Ciphertext bootstrapping benchmark for the revised parameter sets.
//
//   ./bench_boot set=po2-80 mode=on trials=3
//   ./bench_boot list
//
// Reports, per run: the chain HElib actually built, the slot count, the time
// of each bootstrapping stage separately, the capacity remaining, and whether
// decryption round-tripped. The stage breakdown is the point: the contribution
// changes digit extraction only, so a total-time ratio understates it and a
// digit-extraction ratio alone overstates what a user sees. Both are printed.
//
// mode=off selects the bounded-support digit-extraction evaluator, mode=on the
// order-four evaluator. Everything else in the pipeline is identical, so the
// ratio between two runs isolates the contribution.

#include <chrono>
#include <cstdlib>
#include <cstring>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

#include <helib/helib.h>
#include <helib/matmul.h>
#include <helib/fhe_stats.h>

#include "parameter.h"

using namespace order4;

namespace {

struct Args {
  std::string set = "po2-80";
  std::string mode = "on";
  long trials = 3;
  long p_override = 0;   // pick a different plaintext prime from PO2_PRIMES
  bool list = false;
};

Args parse(int argc, char** argv) {
  Args a;
  for (int i = 1; i < argc; ++i) {
    std::string s(argv[i]);
    if (s == "list") { a.list = true; continue; }
    auto eq = s.find('=');
    if (eq == std::string::npos) continue;
    std::string k = s.substr(0, eq), v = s.substr(eq + 1);
    if (k == "set") a.set = v;
    else if (k == "mode") a.mode = v;
    else if (k == "trials") a.trials = std::stol(v);
    else if (k == "p") a.p_override = std::stol(v);
  }
  return a;
}

const ParamSet* find(const std::string& name) {
  for (size_t i = 0; i < N_PO2; ++i)
    if (name == PO2_SETS[i].name) return &PO2_SETS[i];
  for (size_t i = 0; i < N_GEN; ++i)
    if (name == GEN_SETS[i].name) return &GEN_SETS[i];
  return nullptr;
}

void listSets() {
  std::cout << "set          p       m        phi(m)  log2Q   slots   "
               "lambda_stock  lambda_after  attack in scope\n";
  auto row = [](const ParamSet& s) {
    std::cout << std::left << std::setw(12) << s.name << " "
              << std::right << std::setw(6) << s.p << "  "
              << std::setw(7) << s.m << "  " << std::setw(6) << s.phim << "  "
              << std::setw(6) << (long)s.log2Q << "  " << std::setw(6) << s.slots
              << "  " << std::setw(11) << s.lambda_stock
              << "  " << std::setw(12) << s.lambda_after
              << "  " << (s.power_of_two ? "yes" : "no") << "\n";
  };
  for (size_t i = 0; i < N_PO2; ++i) row(PO2_SETS[i]);
  for (size_t i = 0; i < N_GEN; ++i) row(GEN_SETS[i]);
  std::cout << "\nPlaintext primes measured on the power-of-two 80 bit set, "
               "with their order-four radix:\n  ";
  for (size_t i = 0; i < N_PO2_PRIME; ++i)
    std::cout << PO2_PRIMES[i].p << "(A=" << PO2_PRIMES[i].A << ") ";
  std::cout << "\n\nWithdrawn as recommendations, kept only as timing "
               "demonstrations:\n  ";
  for (size_t i = 0; i < N_WITHDRAWN; ++i)
    std::cout << WITHDRAWN[i].name << "(" << WITHDRAWN[i].lambda_after << ") ";
  std::cout << "\n";
}

double now() {
  using namespace std::chrono;
  return duration<double>(steady_clock::now().time_since_epoch()).count();
}

}  // namespace

int main(int argc, char** argv) {
  Args a = parse(argc, argv);
  if (a.list) { listSets(); return 0; }

  const ParamSet* ps = find(a.set);
  if (!ps) {
    std::cerr << "unknown set '" << a.set << "'; run `bench_boot list`\n";
    return 1;
  }
  if (a.mode != "on" && a.mode != "off") {
    std::cerr << "mode must be on or off\n";
    return 1;
  }

  long p = ps->p, A = ps->A;
  if (a.p_override) {
    bool found = false;
    for (size_t i = 0; i < N_PO2_PRIME; ++i)
      if (PO2_PRIMES[i].p == a.p_override) { p = PO2_PRIMES[i].p; A = PO2_PRIMES[i].A; found = true; }
    if (!found) {
      std::cerr << "p=" << a.p_override << " is not in the measured table. "
                   "Its order-four radix would be unknown, and passing another "
                   "prime's radix silently disables the order-four path, which "
                   "produces a run that measures nothing.\n";
      return 1;
    }
  }

  // The evaluator is selected by environment variables read inside the library.
  // Setting them here keeps the choice next to the parameter set rather than in
  // the caller's shell, where a mismatch is easy to miss.
  setenv("HELIB_AUX_ORDER4_EVAL", a.mode == "on" ? "1" : "0", 1);
  if (a.mode == "on") setenv("HELIB_EXPLICIT_AUX", std::to_string(A).c_str(), 1);

  std::cout << "set=" << ps->name << " p=" << p << " A=" << A
            << " m=" << ps->m << " phi(m)=" << ps->phim
            << " h_main=" << ps->h_main << " h_enc=" << ps->h_enc
            << " mode=" << a.mode << " trials=" << a.trials << "\n"
            << "requesting bits=" << ps->bits_request
            << ", the chain that produced log2Q=" << ps->log2Q
            << " on the reference machine\n" << std::flush;

  NTL::Vec<long> mvec;
  std::vector<long> gens, ords;
  if (ps->power_of_two) {
    // A power-of-two ring is a single prime power, and for p = 1 mod 4 the
    // hypercube must have two dimensions: the generator -1 of order 2 together
    // with 5, whose order in the quotient is phi(m)/(2*d). A one dimensional
    // form compiles and builds a context but then fails inside the linear
    // transform, so this is not a free choice.
    mvec.SetLength(1); mvec[0] = ps->m;
    gens = {ps->m - 1, 5};
    ords = {2, ps->slots / 2};
  } else {
    std::cerr << "General cyclotomic sets need their mvec, gens and ords, which "
                 "depend on m = q1*q2 and on ord_m(p). Generate them with "
                 "security/crt_params.py and paste them here; they are not "
                 "derivable from the table alone.\n";
    return 2;
  }

  double t0 = now();
  helib::Context context = helib::ContextBuilder<helib::BGV>()
                               .m(ps->m)
                               .p(p)
                               .r(1)
                               .gens(gens)
                               .ords(ords)
                               .bits(ps->bits_request)
                               .c(3)
                               .bootstrappable(true)
                               .mvec(mvec)
                               .skHwt(ps->h_main)
                               .build();
  double t_ctx = now() - t0;

  const long builtBits = context.bitSizeOfQ();
  std::cout << "built chain = " << builtBits << " bits"
            << " (table records " << ps->log2Q << ")\n"
            << "slots = " << context.getNSlots()
            << ", context setup = " << std::fixed << std::setprecision(2)
            << t_ctx << " s\n" << std::flush;
  if (std::abs((double)builtBits - ps->log2Q) > 40)
    std::cout << "NOTE: the built chain differs from the table by more than 40 "
                 "bits. The `bits` request is not portable across builds; quote "
                 "the built value, not the request.\n";

  t0 = now();
  helib::SecKey secretKey(context);
  secretKey.GenSecKey();
  helib::addSome1DMatrices(secretKey);
  helib::addFrbMatrices(secretKey);
  secretKey.genRecryptData();
  double t_key = now() - t0;
  const helib::PubKey& publicKey = secretKey;
  std::cout << "key generation = " << t_key << " s\n" << std::flush;

  const helib::EncryptedArray& ea = context.getEA();
  std::vector<long> ptxt(ea.size());
  for (size_t i = 0; i < ptxt.size(); ++i) ptxt[i] = (long)(i % p);

  std::vector<double> times;
  bool allCorrect = true;
  for (long t = 0; t < a.trials; ++t) {
    helib::Ctxt c(publicKey);
    ea.encrypt(c, publicKey, ptxt);
    // drive the ciphertext down so that bootstrapping has work to do
    while (c.bitCapacity() > 60) c.multByConstant(NTL::ZZX(1));

    double b0 = now();
    publicKey.reCrypt(c);
    double dt = now() - b0;
    times.push_back(dt);

    std::vector<long> out;
    ea.decrypt(c, secretKey, out);
    bool ok = (out.size() == ptxt.size());
    for (size_t i = 0; ok && i < out.size(); ++i) ok = (out[i] == ptxt[i]);
    allCorrect = allCorrect && ok;

    std::cout << "  trial " << t << ": recrypt " << dt << " s, capacity after "
              << c.bitCapacity() << " bits, decrypt "
              << (ok ? "correct" : "WRONG") << "\n" << std::flush;
  }

  double mean = 0;
  for (double x : times) mean += x;
  mean /= times.size();
  std::cout << "\nRESULT set=" << ps->name << " p=" << p << " mode=" << a.mode
            << " built_log2Q=" << builtBits
            << " slots=" << context.getNSlots()
            << " mean_recrypt_s=" << mean
            << " correct=" << (allCorrect ? 1 : 0) << "\n";
  std::cout << "Per stage timings, if the library was built with timing enabled, "
               "are printed below by helib::printAllTimers.\n";
  helib::printAllTimers(std::cout);
  return allCorrect ? 0 : 3;
}
