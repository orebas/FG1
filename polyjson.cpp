#include "polyjson.hpp"
#include "arbxx.hpp"
#include "benchmark.hpp"
#include <iostream>

struct polyjson {
  std::vector<std::pair<mpfr::mpreal, mpfr::mpreal>> coeffs;
  bool sparse;
  std::vector<long int> sparse_indices;
  std::string filename;
};

void parsePol2(const std::string &polfilename) {
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

int main(int argc, char **argv) {

  for (auto const &dir_entry : std::filesystem::directory_iterator{"."}) {
    if (dir_entry.path().extension() == ".pol") {
      // std::cout << dir_entry << " " << dir_entry.path().extension() <<
      // std::endl;
      parsePol2(dir_entry.path().string());
    }
    // std::cout << dir_entry << " " << dir_entry.path().extension() <<
    // std::endl;
  }
}

