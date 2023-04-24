#pragma once

#include "mpreal.h"
#include <array>
#include <cmath>
#include <complex>
#include <fstream>
#include <iostream>
#include <iterator>
#include <vector>

#include "benchmark.hpp"

#include <nlohmann/json.hpp>
using json = nlohmann::json;

std::pair<std::vector<mpfr::mpreal>, std::vector<mpfr::mpreal>>
PolToMPRealVectors(mps_context *s, mps_monomial_poly *p, slong wp);

namespace mpfr {
void to_json(json &j, const mpreal &p) { j = json(p.toString()); }

void from_json(const json &j, mpreal &p) {
  std::string d;
  j.get_to(d);
  p = mpfr::mpreal(d);
}
}; // namespace mpfr

struct polyjson {
  std::vector<mpfr::mpreal> RealCoefficients;
  std::vector<mpfr::mpreal> ImagCoefficients;

  bool sparse;
  std::vector<long int> sparse_indices;
  std::string filename;
  NLOHMANN_DEFINE_TYPE_INTRUSIVE(polyjson, RealCoefficients, ImagCoefficients,
                                 sparse, sparse_indices, filename);
};

ComplexPoly 
// below returns the name of the file saved, or any empty string on failure.
std::string saveJSON(mps_context *local_s, mps_polynomial *poly_local_poly,
                     const std::string &polfilename) {
  const int desired_prec = 10000;
  mps_monomial_poly *local_poly =
      MPS_POLYNOMIAL_CAST(mps_monomial_poly, poly_local_poly);
  mps_context_set_output_prec(local_s, desired_prec);
  mps_context_set_output_goal(local_s, MPS_OUTPUT_GOAL_APPROXIMATE);
  slong prec = 2000; // TODO(orebas) MAGIC NUMBER
  ComplexPoly acb_style_poly = ComplexPoly(local_s, local_poly, prec);

  auto poly_vectors = PolToMPRealVectors(local_s, local_poly, desired_prec);
  std::cout << poly_vectors.first << poly_vectors.second << std::endl;
  polyjson to_save;
  to_save.RealCoefficients = poly_vectors.first;
  to_save.ImagCoefficients = poly_vectors.second;
  to_save.sparse = false;
  to_save.filename = polfilename + ".json";
  json j1(to_save);
  std::cout << j1 << std::endl;
  std::ofstream fileoutput(to_save.filename);
  fileoutput << j1 << std::endl;
  fileoutput.close();
  return to_save.filename;
}

std::string PolfileToJson(const std::string &polfilename) {
  mps_context *local_s = nullptr;
  mps_polynomial *local_poly = nullptr;
  local_s = mps_context_new();
  slong input_precision = -1;
  char *inline_poly = nullptr;

  FILE *infile = nullptr;
  mps_phase phase = no_phase;

  mps_boolean explicit_algorithm_selection = false;
  infile = fopen(polfilename.c_str(), "re");
  if (infile == nullptr) {
    mps_error(local_s, "Cannot open input file for read, aborting.");
    mps_print_errors(local_s);
    return ""; // EXIT_FAILURE;
  }

  /* Parse the input stream and if a polynomial is given as output,
   * allocate also a secular equation to be used in regeneration */
  local_poly = mps_parse_stream(local_s, infile);

  if (local_poly == nullptr) {
    mps_error(local_s, "Error while parsing the polynomial, aborting.");
    mps_print_errors(local_s);
    return ""; // EXIT_FAILURE;
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

  return saveJSON(local_s, local_poly, polfilename);
}

std::pair<std::vector<mpfr::mpreal>, std::vector<mpfr::mpreal>>
PolToMPRealVectors(mps_context *s, mps_monomial_poly *p, slong wp) {
  mpfr::mpreal::set_default_prec(wp);
  rdpe_t u;
  // cdpe_t cx;
  std::vector<mpfr::mpreal> resultsReal;
  std::vector<mpfr::mpreal> resultsComplex;
  pthread_mutex_lock(&p->mfpc_mutex[0]);
  if (mpc_get_prec(p->mfpc[0]) < wp) {
    std::cout << "it is " << mpc_get_prec(p->mfpc[0]) << " and " << wp
              << std::endl;
    pthread_mutex_unlock(&p->mfpc_mutex[0]);
    mps_monomial_poly_raise_precision(s, MPS_POLYNOMIAL(p), wp);
    std::cout << "it is " << mpc_get_prec(p->mfpc[0]) << " and " << wp
              << std::endl;
  } else {
    pthread_mutex_unlock(&p->mfpc_mutex[0]);
  }

  /* Set 4 * machine precision in u */
  rdpe_set_2dl(u, 1.0, 2 - wp);

  slong local_prec =
      std::max(mpf_get_prec(p->mfpc[0]->r), mpf_get_prec(p->mfpc[0]->i));
  local_prec = std::max(local_prec, wp);
  mpfr_set_default_prec(local_prec);
  mpfr_t mpfrre;
  mpfr_t mpfrim;
  mpfr_init2(mpfrre, local_prec);
  mpfr_init2(mpfrim, local_prec);

  for (int j = 0; j <= MPS_POLYNOMIAL(p)->degree; j++) {
    mpfr_set_f(mpfrre, (p->mfpc[j])->r, MPFR_RNDN);
    mpfr_set_f(mpfrim, (p->mfpc[j])->i, MPFR_RNDN);
    mpfr::mpreal mprealre(mpfrre); // make sure the default precision is ok?
    mpfr::mpreal mprealim(mpfrim);
    std::complex<mpfr::mpreal> coeff;
    resultsReal.push_back(mprealre);
    resultsComplex.push_back(mprealim);
  }
  mpfr_clear(mpfrre);
  mpfr_clear(mpfrim);
  return std::pair<std::vector<mpfr::mpreal>, std::vector<mpfr::mpreal>>(
      resultsReal, resultsComplex);
}