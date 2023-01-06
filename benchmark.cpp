

#include <chrono>
#include <cmath>
#include <iostream>

#include <filesystem>
#include <functional>
#include <memory>
#include <utility>

//#include "mpreal.h" this needs to be incuded before all the arb/flint stuff
#include "acb_calc.h"
#include "acb_poly.h"
#include "arb.h"

#include "acb.h"
#include "arb_fmpz_poly.h"
#include "flint/arith.h"
#include "flint/profiler.h"
#include "tabulate.hxx"

#define _MPS_PRIVATE
#include <cstring>
#include <mps/mps.h>

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

constexpr auto MPSOLVE_GETOPT_STRING = "a:G:D:d::t:o:O:j:S:O:i:vl:bp:rs:c";

#ifndef __WINDOWS
#include <csignal>

void status(int signal, mps_context *mps_c);

#undef _MPS_PRIVATE
#endif

#include <benchmark.hpp>
// #include "Hungarian.h"
//  void * cleanup_context(mps_context *ctx, void *user_data);

void timeSolvers(std::vector<ComplexPoly> pvec);

int main(int argc, char **argv) {
  // mpfr_set_default_prec(3000);
  //  swap debugging
  /*
      ACB a(12, 34, 120),
          b(56, 78, 240);

      std::cout << a << " " << b;
      std::cout << a.intprec << " " << b.intprec << std::endl;
      // std::swap(a, b);

      ACB temp = std::move(a);  // or T temp(std::move(t1));

      std::cout << temp << " " << b << " ";
      std::cout << " " << temp.intprec << " " << b.intprec << " " << std::endl;

      a = std::move(b);

      std::cout << temp << " " << a << " ";
      std::cout << " " << temp.intprec << " " << a.intprec << " " << std::endl;

      b = std::move(temp);

      std::cout << a << " " << b << " ";
      std::cout << " " << a.intprec << " " << b.intprec << " " << std::endl;

      // std::cout << a << " " << b << " " << temp;
      // std::cout << " " << a.intprec << " " << b.intprec << " " <<
     temp.intprec << std::endl;

      // std::cout << a << " " << b << " " << temp;
      // std::cout << " " << a.intprec << " " << b.intprec << " " <<
     temp.intprec << std::endl;

      //    std::cout << " " << a << " " << b;
      //   std::cout << " " << a.intprec << " " << b.intprec << " " <<
     std::endl;
  */
  // end swap debugging
  std::vector<ComplexPoly> pvec;
  std::vector<int> dvec;
  const int mindegree = 5;
  const int maxdegree = 40;
  const int inputprecision = 512;
  for (int i = 1; i <= maxdegree; i++) {
    dvec.push_back(i);
    if (i >= mindegree) {
      pvec.push_back(polyFromRoots<int>(
          dvec, inputprecision)); // TODO(orebas) MAGIC NUMBER
    }
  }

  timeSolvers(pvec);
  // auto res = p.MPSolve(120);
  return 0;
  for (auto const &dir_entry : std::filesystem::directory_iterator{"."}) {
    if (dir_entry.path().extension() == ".pol") {
      // std::cout << dir_entry << " " << dir_entry.path().extension() <<
      // std::endl;
      parsePol(dir_entry.path().string());
    }
    // std::cout << dir_entry << " " << dir_entry.path().extension() <<
    // std::endl;
  }
}

class Solver { // TODO(orebas) add RAII make a real class
public:
  std::string name;
  std::function<std::vector<ACB>(ComplexPoly &)> solver;
};

void timeSolvers(std::vector<ComplexPoly> pvec) {
  SimpleNewtonSolver sns;
  const slong prec = 512; // TODO(orebas) make this a parameter, for each poly.
  Solver mpsolve;
  mpsolve.name = "MPSolve";
  mpsolve.solver = [](const ComplexPoly &poly) -> std::vector<ACB> {
    return poly.MPSolve(prec);
  };

  Solver arblibsolve;
  arblibsolve.name = "ArbLib";
  arblibsolve.solver = [](const ComplexPoly &poly) -> std::vector<ACB> {
    return poly.ArbLibSolve(prec);
  };
  Solver simplenewton;
  simplenewton.name = "Simple Newton";
  simplenewton.solver =
      [&sns](ComplexPoly &poly) -> std::vector<ACB> { // try to const this later
    return sns.solve(poly, prec);
  };
  std::vector<Solver> solvers;
  solvers.push_back(mpsolve);
  solvers.push_back(arblibsolve);
  solvers.push_back(simplenewton);
  Array3d<double> timings(pvec.size(), solvers.size(), 1);
  Array3d<std::vector<ACB>> results(pvec.size(), solvers.size(), 1,
                                    std::vector<ACB>());
  // std::vector<ACB> results[pvec.size()][solvers.size()];
  for (std::size_t i = 0; i < pvec.size(); i++) {
    for (std::size_t j = 0; j < solvers.size(); j++) {
      results(i, j, 0) = solvers[j].solver(pvec[i]);
      timings(i, j, 0) = measure<std::chrono::milliseconds>::execution(
          solvers[j].solver, pvec[i]);
    }
  }
  Array3d<ARB> resultsdiff(pvec.size(), solvers.size(), 1, ARB(0, prec));
  for (std::size_t i = 0; i < pvec.size(); i++) {
    for (std::size_t j = 0; j < solvers.size(); j++) {
      resultsdiff(i, j, 0) = L2Norm(results(i, j, 0), results(i, 0, 0), prec);
    }
  }
  using namespace tabulate;
  Table output_table;
  const int digits = 20;
  for (std::size_t i = 0; i < pvec.size(); i++) {
    output_table.add_row({
        std::to_string(i),
        solvers[0].name,
        std::to_string(timings(i, 0, 0)),
        arb_get_str(resultsdiff(i, 0, 0).r, digits, ARB_STR_NO_RADIUS),
        solvers[1].name,
        std::to_string(timings(i, 1, 0)),
        arb_get_str(resultsdiff(i, 1, 0).r, digits, ARB_STR_NO_RADIUS),
        solvers[2].name,
        std::to_string(timings(i, 2, 0)),
        arb_get_str(resultsdiff(i, 2, 0).r, digits, ARB_STR_NO_RADIUS),

    });
  }
  std::cout << output_table << std::endl;
}

void *cleanup_context(mpscp ctx, void *user_data, mps_polynomial *poly,
                      mpscp &s_to_zero, bool outp) {
  /* Check for errors */
  if (mps_context_has_errors(ctx)) {
    mps_print_errors(ctx);
    return nullptr;
  }

  /* Output the roots */
  if (!outp) {
    mps_output(ctx);
  }
  /* Free used data */
  if (poly != nullptr) {
    mps_polynomial_free(ctx, poly);
  }
  mps_context_free(ctx);

  s_to_zero = nullptr;

  return nullptr;
}

void PolfileToJson(const std::string &polfilename) {
  mps_context *local_s = nullptr;
  mps_polynomial *local_poly = nullptr;
  slong input_precision = -1;
  char *inline_poly = nullptr;

  FILE *infile = nullptr;
  mps_phase phase = no_phase;

  mps_boolean explicit_algorithm_selection = false;
  infile = fopen(polfilename.c_str(), "re");
  if (infile == nullptr) {
    mps_error(local_s, "Cannot open input file for read, aborting.");
    mps_print_errors(local_s);
    return; // EXIT_FAILURE;
  }

  /* Parse the input stream and if a polynomial is given as output,
   * allocate also a secular equation to be used in regeneration */
  local_poly = mps_parse_stream(local_s, infile);

  if (local_poly == nullptr) {
    mps_error(local_s, "Error while parsing the polynomial, aborting.");
    mps_print_errors(local_s);
    return; // EXIT_FAILURE;
  } else {
    mps_context_set_input_poly(local_s, local_poly);
  }

  if (input_precision >= 0) {
    mps_polynomial_set_input_prec(local_s, local_poly,
                                  2000); // TODO(orebas) MAGIC NUMBER
  }
  if (!explicit_algorithm_selection) {
    mps_context_select_algorithm(local_s,
                                 (MPS_IS_MONOMIAL_POLY(local_poly) &&
                                  MPS_DENSITY_IS_SPARSE(local_poly->density))
                                     ? MPS_ALGORITHM_STANDARD_MPSOLVE
                                     : MPS_ALGORITHM_SECULAR_GA);
  }
  /* Close the file if it's not stdin */
  if (inline_poly == nullptr) {
    (void)fclose(infile);
  }

  mps_context_set_starting_phase(local_s, phase);

  saveJSON(local_s, local_poly, polfilename);
}

void parsePol(const std::string &polfilename) {
  mps_context *local_s = nullptr;
  mps_polynomial *local_poly = nullptr;

  /* Create a new status */
  local_s = mps_context_new();

  /* Associate the SIGUSR1 signal to the mps_dump () function */
  // #ifndef __WINDOWS
  //     signal(SIGUSR1, status);
  // #endif

  /* Leave this value to -1 if not precision has been enforced
   * by the user. Otherwise, override the polynomial input
   * precision. */
  slong input_precision = -1;
  char *inline_poly = nullptr;

  FILE *infile = nullptr;

  /* Parse options */
  // mps_opt *opt = nullptr; //commented out but we might need it later
  mps_phase phase = no_phase;

  /* This will be true if the user has explicitely selected the algorithm,
     otherwise we will be using our own heuristic. */
  mps_boolean explicit_algorithm_selection = false;
  // opt = nullptr;
  /* we could set a logfile with the following code:
 {
              FILE *logstr = fopen(opt->optvalue, "w");
              if (!logstr)
                  mps_error(s, "Cannot open selected log file.");
              mps_context_set_log_stream(s, logstr);
          }
   */
  // we could set output format with mps_context_set_output_format //see
  // main() for examples
  infile = fopen(polfilename.c_str(), "re");
  if (infile == nullptr) {
    mps_error(local_s, "Cannot open input file for read, aborting.");
    mps_print_errors(local_s);
    return; // EXIT_FAILURE;
  }

  /* Parse the input stream and if a polynomial is given as output,
   * allocate also a secular equation to be used in regeneration */
  local_poly = mps_parse_stream(local_s, infile);

  if (local_poly == nullptr) {
    mps_error(local_s, "Error while parsing the polynomial, aborting.");
    mps_print_errors(local_s);
    return; // EXIT_FAILURE;
  } else {
    mps_context_set_input_poly(local_s, local_poly);
  }
  /*local_s->active_poly = local_poly;
  std::cout << local_s << std::endl;
  std::cout << local_s->active_poly << std::endl;
  std::cout << local_poly << std::endl;
*/
  /* Override input precision if needed */
  if (input_precision >= 0) {
    mps_polynomial_set_input_prec(local_s, local_poly,
                                  2000); // TODO(orebas) MAGIC NUMBER
  }
  /* Perform some heuristic for the algorithm selection, but only if the
   * user hasn't explicitely selected one. */
  if (!explicit_algorithm_selection) {
    mps_context_select_algorithm(local_s,
                                 (MPS_IS_MONOMIAL_POLY(local_poly) &&
                                  MPS_DENSITY_IS_SPARSE(local_poly->density))
                                     ? MPS_ALGORITHM_STANDARD_MPSOLVE
                                     : MPS_ALGORITHM_SECULAR_GA);
  }

  /* Close the file if it's not stdin */
  if (inline_poly == nullptr) {
    (void)fclose(infile);
  }

  /* Select the starting phase according to user input */
  mps_context_set_starting_phase(local_s, phase);

  // BEGIN WILK20 DEBUGGING CODE
  /*
  mps_monomial_poly *local_poly2 = MPS_POLYNOMIAL_CAST(mps_monomial_poly,
  local_poly);

  fmpz_poly_t f, g;
  fmpz_poly_factor_t fac;
  fmpz_t t;
  acb_ptr roots;
  slong compd, printd, i, j, deg;
  int flags;
  flags = ARB_FMPZ_POLY_ROOTS_VERBOSE;

  fmpz_poly_init(f);
  fmpz_poly_init(g);
  fmpz_init(t);
  fmpz_poly_one(f);

  slong n = 20;
  fmpz_poly_zero(g);
  fmpz_poly_fit_length(g, n + 2);
  arith_stirling_number_1_vec(g->coeffs, n + 1, n + 2);
  _fmpz_poly_set_length(g, n + 2);
  fmpz_poly_shift_right(g, g, 1);
  fmpz_poly_mul(f, f, g);
  i++;

  fmpz_poly_factor_init(fac);
  std::cout << std::endl;
  std::cout << "FMPZ_POLY" << std::endl;
  fmpz_poly_print(f);
  std::cout << std::endl;
  fmpz_poly_fprint_pretty(stdout, f, "x");
  std::cout << std::endl;
  flint_printf("computing squarefree factorization...\n");
  slong prec = 6000;

  ComplexPoly testp(prec), diffp(prec);
  ComplexPoly acb_style_poly = ComplexPoly(local_s, local_poly2, prec);
  acb_poly_set_fmpz_poly(testp.pol, f, prec);
  acb_poly_sub(diffp.pol, testp.pol, acb_style_poly.pol, prec);

  std::cout << "TESTP" << std::endl;
  testp.print(30);

  std::cout << "ACB_STYLE_POLY" << std::endl;
  acb_style_poly.print(30);

  std::cout << "DIFF" << std::endl;
  diffp.print(30);

  TIMEIT_ONCE_START
  fmpz_poly_factor_squarefree(fac, f);
  TIMEIT_ONCE_STOP

  TIMEIT_ONCE_START
  for (i = 0; i < fac->num; i++) {
      deg = fmpz_poly_degree(fac->p + i);

      flint_printf("%wd roots with multiplicity %wd\n", deg, fac->exp[i]);
      roots = _acb_vec_init(deg);

      arb_fmpz_poly_complex_roots(roots, fac->p + i, flags, compd * 3.32193
  + 2);

      if (printd) {
          for (j = 0; j < deg; j++) {
              acb_printn(roots + j, printd, 0);
              flint_printf("\n");
          }
      }

      _acb_vec_clear(roots, deg);
  }
  TIMEIT_ONCE_STOP

  fmpz_poly_factor_clear(fac);
  fmpz_poly_clear(f);
  fmpz_poly_clear(g);
  fmpz_clear(t);

  // END WILK20 DEBUGGING CODE
  */
  solveCompare(local_s, local_poly, polfilename);
}

void temp_print(mps_context *local_s, mps_polynomial *poly) {
  // mps_monomial_poly *p = MPS_MONOMIAL_POLY(poly);
  mps_monomial_poly *mon_poly = reinterpret_cast<mps_monomial_poly *>(poly);

  std::cout << poly->degree << std::endl;
  (void)fflush(stdout);
  const int digits = 10;
  for (int i = 0; i <= poly->degree; i++) {
    rdpe_out(mon_poly->dap[i]);
    std::cout << mon_poly->fap[i] << std::endl;
    std::cout << mon_poly->fpr[i] << std::endl;
    cplx_out(mon_poly->fpc[i]);
    rdpe_out(mon_poly->dpr[i]);
    cdpe_out(mon_poly->dpc[i]);
    mpf_out_str(stdout, digits, 0, mon_poly->mfpr[i]);
    mpc_out_str(stdout, digits, 0, mon_poly->mfpc[i]);

    if (i < poly->degree) {
      mpc_out_str(stdout, digits, 0, mon_poly->mfppc[i]);
    }
    mpq_out_str(stdout, digits, mon_poly->initial_mqp_r[i]);
    mpq_out_str(stdout, digits, mon_poly->initial_mqp_i[i]);
    std::cout << mon_poly->spar[i] << std::endl;
  }
}

std::vector<ACB> powersums(slong depth, slong max_k, slong prec) {
  std::vector<ACB> results;
  for (int d = 0; d < depth + 1; d++) {
    ACB temp(0, 0, prec);

    for (int i = 1; i <= max_k; i++) {
      temp += (ACB(i, 0, prec)).power(d);
    }
    results.push_back(temp);
  }
  return results;
}

void saveJSON(mps_context *local_s, mps_polynomial *poly_local_poly,
              const std::string &polfilename) {
  mps_monomial_poly *local_poly =
      MPS_POLYNOMIAL_CAST(mps_monomial_poly, poly_local_poly);
  mps_context_set_output_prec(local_s, desired_prec);
  mps_context_set_output_goal(local_s, MPS_OUTPUT_GOAL_APPROXIMATE);
  slong prec = 2000; // TODO(orebas) MAGIC NUMBER
  ComplexPoly acb_style_poly = ComplexPoly(local_s, local_poly, prec);

}

void solveCompare(mps_context *local_s, mps_polynomial *poly_local_poly,
                  const std::string &polfilename) {
  mps_monomial_poly *local_poly =
      MPS_POLYNOMIAL_CAST(mps_monomial_poly, poly_local_poly);

  std::cout << "Starting " << polfilename << std::endl;
  const int desired_prec = 512;
  mps_context_set_output_prec(local_s, desired_prec);
  mps_context_set_output_goal(local_s, MPS_OUTPUT_GOAL_APPROXIMATE);
  /* Solve the polynomial */
  auto funcMPS = [&]() -> void { mps_mpsolve(local_s); };

  slong prec = 2000; // TODO(orebas) MAGIC NUMBER
  ComplexPoly acb_style_poly = ComplexPoly(local_s, local_poly, prec);
  std::cout << std::endl;
  std::vector<ComplexPoly> pvec;
  std::vector<ComplexPoly> qvec;
  // DLGFG(acb_style_poly, pvec, qvec, 90);
  const int depth = 20;
  // FGSolver FGS(acb_style_poly, depth);
  // FGS.generatePnQnPolynomials();
  /*std::cout << "PVEC:";
  for (const auto &p : FGS.pvec) {
      p.print();
  }
  std::cout << "QVEC:";

  for (const auto &q : FGS.qvec) {
      q.print();
  }*/

  // std::vector<std::vector<ARB>> radii(depth);
  // std::vector<std::vector<ACB>> rootapprox(depth);

  // FGS.DLGRadii(radii);
  //  FGS.FGRoots(rootapprox);

  std::vector<ACB> RootSolver1Roots;
  RootSolver1 RS1(acb_style_poly, depth);
  auto funcRS1 = [&]() -> void { RootSolver1Roots = RS1.solve(); };
  // std::cout << radii;
  /*std::cout << "Root radii array is ";
  std::cout << radii;

  std::cout << "Root approx array is ";
  std::cout << rootapprox;*/
  //  std::cout << powersums(110, 20, 6000);
  // std::cout    << "Begin RRE debugging" << std::endl;
  // RecursiveRadiiEstimator RRE(acb_style_poly);

  /*for (int d = 0; d < 15; d++) {
      std::cout << "RRE " << d << " " << std::endl;
      FGS.pvec[d].print();
      for (int n = 0; n <= acb_style_poly.degree(); n++) {
          std::cout << RRE.calcCoeff(d, n);
      }
  }*/

  /*for (int depth = 0; depth < 15; depth++) {
      std::cout << depth << "[min]" << RRE.minRadius(depth);
      std::cout << "[max]" << RRE.maxRadius(depth);
  }

  for (int depth = 0; depth < 15; depth++) {
      std::cout << depth << "[minRoot]" << RRE.minRoot(depth);

      std::cout << depth << "[minRoot]" << RRE.maxRoot(depth);
  }*/
  // uses recursive FG and DLG to estimate largest and smallest roots, and
  // root radii

  double mpsTime =
      measure<std::chrono::milliseconds>::execution(funcMPS); // runs mpsolve
  double RS1Time =
      measure<std::chrono::milliseconds>::execution(funcRS1); // runs mpsolve

  (void)fflush(stdout);
  mpc_t *results = nullptr;
  // std::cout << "MPS SOLVE" << std::endl;
  mps_context_get_roots_m(local_s, &results, nullptr);
  std::vector<ACB> mps_roots;
  for (slong rt = 0; rt < acb_style_poly.degree(); rt++) {
    mps_roots.emplace_back(ACB(results[rt], prec));
  }

  // std::cout << "MPS SOLVE1" << std::endl;
  //  std::cout << mps_roots;
  struct {
    bool operator()(const ACB &a, const ACB &b) const {
      // ARF l = a.abs_ubound_arf();
      // ARF r = b.abs_ubound_arf();
      // int i = arf_cmp(l.f, r.f);
      // return i < 0;
      if (arb_lt(a.abs().r, b.abs().r) != 0) {
        return true;
      }
      if (arb_gt(a.abs().r, b.abs().r) != 0) {
        return false;
      }
      ARB arg1(0.0, a.intprec);
      ARB arg2(0.0, b.intprec);
      acb_arg(arg1.r, a.c, a.intprec);
      acb_arg(arg2.r, b.c, a.intprec);
      if (arb_lt(arg1.r, arg2.r) != 0) {
        return true;
      }
      if (arb_gt(arg1.r, arg2.r) != 0) {
        return false;
      }
      return false;
    }
  } customLess;

  struct {
    bool operator()(const ACB &a, const ACB &b) const {
      // ARF l = a.abs_ubound_arf();
      // ARF r = b.abs_ubound_arf();
      // int i = arf_cmp(l.f, r.f);
      // return i < 0;

      ARB re1(0.0, a.intprec);
      ARB im1(0.0, b.intprec);

      ARB re2(0.0, a.intprec);
      ARB im2(0.0, b.intprec);

      acb_get_real(re1.r, a.c);
      acb_get_real(re2.r, b.c);
      acb_get_imag(im1.r, a.c);
      acb_get_imag(im2.r, b.c);

      if (arb_lt(re1.r, re2.r)) {
        return true;
      }
      if (arb_gt(re1.r, re2.r)) {
        return false;
      }
      if (arb_lt(im1.r, im2.r)) {
        return true;
      }
      if (arb_gt(im1.r, im2.r)) {
        return false;
      }
      return false;
    }
  } customLessLex;

  std::sort(mps_roots.begin(), mps_roots.end(), customLessLex);
  std::sort(RootSolver1Roots.begin(), RootSolver1Roots.end(), customLess);

  std::vector<std::vector<ARB>> costMatrix;
  slong matsize = mps_roots.size();
  costMatrix.resize(matsize);
  // HungarianAlgorithm<ARB> HungAlgo;
  /*  COMMENTED OUT HUNGARIAN ALGORITHM
  for (int i = 0; i < matsize; i++) {
      for (int j = 0; j < matsize; j++) {
          costMatrix[i].push_back((mps_roots[i] -
  RootSolver1Roots[j]).abs());
      }
  }
  std::vector<int> assignment;
  ARB epsilon(0, prec);
  bool epsset = false;

  for (int i = 0; i < matsize; i++) {
      for (int j = 0; j < matsize; j++) {
          if (costMatrix[i][j] > ARB(0, 6000)) {
              if (epsset) {
                  epsilon = costMatrix[i][j];
              } else {
                  if (costMatrix[i][j] < epsilon) {
                      epsilon = costMatrix[i][j];
                  }
              }
          }
      }
  }
  epsilon = MIN(epsilon, ARB(DBL_EPSILON, 6000));
  epsilon = epsilon * epsilon;

  ARB cost = HungAlgo.Solve(costMatrix, assignment, epsilon, ARB(0, prec));

  std::cout << assignment;*/

  std::cout << "RootSolver1" << std::endl;
  std::cout << RootSolver1Roots;

  std::cout << "MPS SOLVE2" << std::endl;
  std::cout << mps_roots;

  std::cout << polfilename << " took " << mpsTime << "ms. for MPS "
            << std::endl;
  std::cout << polfilename << " took " << RS1Time << "ms. for RS1."
            << std::endl;
  // auto res = ACBVectorComp(mps_roots, rootapprox.back());
  // std::cout << res << std::endl;
  std::cout << mps_roots.size() << std::endl;

  auto rootapprox = RS1.rootapprox;

  std::cout << rootapprox.back().size() << std::endl;

  for (std::size_t ri = 1; ri < rootapprox.size(); ri++) {
    std::sort(rootapprox[ri].begin(), rootapprox[ri].end(), customLessLex);
    std::cout << "\nRow " << ri << " "
              << ACBVectorComp(mps_roots, rootapprox[ri]);
  }
  mpc_vfree(results);
  cleanup_context(local_s, nullptr, MPS_POLYNOMIAL(local_poly), local_s, false);
}

void solveDebug(mps_context *local_s, mps_polynomial *poly_local_poly,
                const std::string &polfilename) {
  mps_monomial_poly *local_poly =
      MPS_POLYNOMIAL_CAST(mps_monomial_poly, poly_local_poly);
  /* Solve the polynomial */
  auto func = [&]() -> void { mps_mpsolve(local_s); };

  // ComplexPoly new_poly;
  //  auto parsed_poly = ComplexPoly::from_mps_poly(mps_polynomial *
  //  local_poly);
  slong prec = 2000; // TODO(orebas) MAGIC NUMBER
  ComplexPoly acb_style_poly = ComplexPoly(local_s, local_poly, prec);
  std::vector<ComplexPoly> pvec, qvec;
  // DLGFG(acb_style_poly, pvec, qvec, 90);
  int depth = 40;
  FGSolver FGS(acb_style_poly, depth);
  FGS.generatePnQnPolynomials();
  std::vector<std::vector<ARB>> radii(depth);
  std::vector<std::vector<ACB>> rootapprox(depth);
  FGS.DLGRadii(radii);
  FGS.FGRoots(rootapprox);
  std::cout << radii;
  std::cout << rootapprox;
  // std::cout << powersums(110, 20, 6000);
  std::cout << "Begin RRE debugging" << std::endl;
  RecursiveRadiiEstimator RRE(acb_style_poly);

  /*for (int d = 0; d < 15; d++) {
      std::cout << "RRE " << d << " " << std::endl;
      FGS.pvec[d].print();
      for (int n = 0; n <= acb_style_poly.degree(); n++) {
          std::cout << RRE.calcCoeff(d, n);
      }
  }*/

  for (int depth = 0; depth < 15; depth++) {
    std::cout << depth << "[min]" << RRE.minRadius(depth);
    std::cout << "[max]" << RRE.maxRadius(depth);
  }

  for (int depth = 0; depth < 15; depth++) {
    std::cout << depth << "[minRoot]" << RRE.minRoot(depth);

    std::cout << depth << "[minRoot]" << RRE.maxRoot(depth);
  }

  double initTime = measure<std::chrono::milliseconds>::execution(func);
  (void)fflush(stdout);
  std::cout << polfilename << " took " << initTime << "ms." << std::endl;
  cleanup_context(local_s, nullptr, MPS_POLYNOMIAL(local_poly), local_s, false);
}

class integrand_parameters {
public:
  ComplexPoly p;
  ARB r;
  slong xpower;
  integrand_parameters(const ComplexPoly &initpol, ARB radius, slong pow)
      : p(initpol.intprec), r(0, radius.intprec), xpower(pow) {
    p.explicitCopy(initpol);
    arb_set(r.r, radius.r);
  }
};

int pprime_over_p_integrand(acb_ptr out, const acb_t inp, void *paramsv,
                            slong order, slong prec) {
  // this is the integrand to evaluate integral p'/p, order-th derivative.
  // inp goes from 0 to 1.

  ACB a(0, 0, prec), b(0, 0, prec), z(0, 0, prec);
  ACB inp_copy(0, 0, prec);
  integrand_parameters *params =
      reinterpret_cast<integrand_parameters *>(paramsv);
  acb_set(inp_copy.c, inp);
  acb_mul_2exp_si(inp_copy.c, inp_copy.c, 1); // inp *=2;
  acb_exp_pi_i(z.c, inp_copy.c,
               inp_copy.intprec);                // z = exp(2 * pi *i * inp)
  acb_mul_arb(z.c, z.c, params->r.r, z.intprec); // z = r e ^(2 * pi * inp);

  acb_poly_evaluate2(a.c, b.c, params->p.pol, z.c, params->p.intprec);
  acb_div(out, b.c, a.c, b.intprec); // out = p'/p at re^(2 pi i inp)
  acb_pow_si(z.c, z.c, params->xpower, z.intprec);
  acb_div(out, out, z.c, z.intprec); // return p'(z)/p(z)z^n

  return 0;
}

// below approximates first depth derivatives of pprime/p via cauchy
// integral formula
void DerivEst(std::vector<ACB> &results, const ComplexPoly &p, ARB minrad,
              std::size_t depth) {
  mag_t tol;
  acb_calc_integrate_opt_t options;

  acb_calc_integrate_opt_init(options);

  slong prec = p.intprec; // p.intprec;
  slong goal = prec;
  mag_init(tol);
  mag_set_ui_2exp_si(tol, 1, -prec);

  ACB deriv_est(0, 0, prec);
  results.resize(depth + 1, ACB(0, 0, prec));
  integrand_parameters params(p, minrad, 0);
  ACB start(0, 0, prec);
  ACB end(1, 0, prec);
  auto func(&pprime_over_p_integrand);
  for (std::size_t k = 0; k <= depth; k++) {
    params.xpower = k;
    acb_calc_integrate(deriv_est.c, func, &params, start.c, end.c, goal, tol,
                       options, prec);
    fmpz_t fact;
    fmpz_fac_ui(fact, k);
    acb_mul_fmpz(deriv_est.c, deriv_est.c, fact, deriv_est.intprec);
    results[k].explicitCopy(deriv_est);
  }
}
