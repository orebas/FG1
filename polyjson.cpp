#include "polyjson.hpp"
#include "arbxx.hpp"
#include "benchmark.hpp"
#include <iostream>

std::vector<ACB> MPSolvePolFile(const std::string &polfilename, slong prec);
std::vector<ACB> solveMPSContext(mps_context *local_s,
                                 mps_polynomial *poly_local_poly,
                                 const std::string &polfilename, slong prec);

void parsePol2(const std::string &polfilename) {
 // std::cout << "test\n";
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
                                  20000); // TODO(orebas) MAGIC NUMBER
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

  mps_monomial_poly *plocal_poly =
      MPS_POLYNOMIAL_CAST(mps_monomial_poly, local_poly);

  slong prec = 2000; // TODO(orebas) MAGIC NUMBER
  ComplexPoly acb_style_poly = ComplexPoly(local_s, plocal_poly, prec);
}

void deleteme() {
  std::random_device rd; // TODO(orebas) make this a static object and only
                         // have one generator
  std::mt19937 mt(rd());
  std::uniform_real_distribution<double> dist(-1.0, 1.0);

  const int d = 10;
  const int tests = 30;
  const int prec = 200;
  std::vector<std::vector<ACB>> bmatrix(tests,
                                        std::vector<ACB>(0, ACB(0, 0, prec)));
  for (int i = 0; i < tests; i++) {
    for (int j = 0; j < d; j++) {
      //auto xd = dist(rd);
      ACB cc(dist(rd), dist(rd), prec);
      bmatrix[i].push_back(cc);
    }
    ARB norm = L1Norm(bmatrix[i], prec);
    for (int j = 0; j < d; j++) {
      bmatrix[i][j] /= ACB(norm);
    }
  }
  //std::cout << bmatrix << std::endl;
}

void compareVectors(std:: string filename, std::vector<ACB> r1, std::vector<ACB> r2, slong prec);

int main(int argc, char **argv) {
  // deleteme();
  for (auto const &dir_entry : std::filesystem::directory_iterator{"."}) {
    if (dir_entry.path().extension() == ".pol") {
      slong prec = 4000;
      std::string jsonFileName = PolfileToJson(dir_entry.path().string());
      if (jsonFileName.empty()) {
        std::cout << "Error converting " << dir_entry.path().string()
                  << std::endl;
        break;
      }
      std::vector<ACB> MPSRoots =
          MPSolvePolFile(dir_entry.path().string(), prec);
      /*polyjson p = loadJson(jsonFileName);
      ComplexPoly poly = toComplexPoly(p);
      std::vector<ACB> rootsFromJson = poly.MPSolve(prec);
      std::sort(rootsFromJson.begin(), rootsFromJson.end(),ACB::customLess);
      std::sort(MPSRoots.begin(), MPSRoots.end(), ACB::customLess);
      assert(rootsFromJson.size() == MPSRoots.size());
      
      compareVectors(jsonFileName, MPSRoots, rootsFromJson, prec);
      */
    }
  }
}


void compareVectors(std:: string filename, std::vector<ACB> r1, std::vector<ACB> r2, slong prec){
fmt::print("Solving {}", filename);  

}




std::vector<ACB> MPSolvePolFile(const std::string &polfilename, slong prec) {
  std::vector<ACB> return_empty;

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
    return return_empty; // EXIT_FAILURE;
  }

  /* Parse the input stream and if a polynomial is given as output,
   * allocate also a secular equation to be used in regeneration */
  local_poly = mps_parse_stream(local_s, infile);

  if (local_poly == nullptr) {
    mps_error(local_s, "Error while parsing the polynomial, aborting.");
    mps_print_errors(local_s);
    return return_empty; // EXIT_FAILURE;
  } else {
    mps_context_set_input_poly(local_s, local_poly);
  }
  /* Override input precision if needed */
  if (input_precision >= 0) {
    mps_polynomial_set_input_prec(local_s, local_poly,
                                  prec); // TODO(orebas) MAGIC NUMBER
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
 
  auto retvec  =  solveMPSContext(local_s, local_poly, polfilename, prec);
  return retvec;
}

std::vector<ACB> solveMPSContext(mps_context *local_s,
                                 mps_polynomial *poly_local_poly,
                                 const std::string &polfilename, slong prec) {
  mps_monomial_poly *local_poly =
      MPS_POLYNOMIAL_CAST(mps_monomial_poly, poly_local_poly);

  //std::cout << "Starting " << polfilename << std::endl;
  const int desired_prec = prec;
  mps_context_set_output_prec(local_s, desired_prec);
  mps_context_set_output_goal(local_s, MPS_OUTPUT_GOAL_APPROXIMATE);
  /* Solve the polynomial */
  auto funcMPS = [&]() -> void { mps_mpsolve(local_s); };
  double mpsTime =
      measure<std::chrono::milliseconds>::execution(funcMPS); // runs mpsolve

  (void)fflush(stdout);
  mpc_t *results = nullptr;
  mps_context_get_roots_m(local_s, &results, nullptr);
  std::vector<ACB> mps_roots;
  for (slong rt = 0; rt < poly_local_poly->degree; rt++) {
    mps_roots.emplace_back(ACB(results[rt], prec));
  }


  std::sort(mps_roots.begin(), mps_roots.end(), ACB::customLessLex);
  //std::cout << polfilename << " took " << mpsTime << "ms. for MPS "
   //         << std::endl;
  //std::cout << mps_roots.size() << std::endl;
  mpc_vfree(results);

  cleanup_context(local_s, nullptr, MPS_POLYNOMIAL(local_poly), local_s, true);

  return mps_roots;
}
