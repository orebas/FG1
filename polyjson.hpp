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

namespace mpfr {
void to_json(json &j, const mpreal &p) { j = json(p.toString()); }

void from_json(const json &j, mpreal &p) {
  std::string d;
  j.get_to(d);
  p = mpfr::mpreal(d);
}
}; // namespace mpfr

void saveJSON(mps_context *local_s, mps_polynomial *poly_local_poly,
              const std::string &polfilename) {
  const int desired_prec = 10000;
  mps_monomial_poly *local_poly =
      MPS_POLYNOMIAL_CAST(mps_monomial_poly, poly_local_poly);
  mps_context_set_output_prec(local_s, desired_prec);
  mps_context_set_output_goal(local_s, MPS_OUTPUT_GOAL_APPROXIMATE);
  slong prec = 2000; // TODO(orebas) MAGIC NUMBER
  ComplexPoly acb_style_poly = ComplexPoly(local_s, local_poly, prec);
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

std::vector<std::complex<mpfr::mpreal>>
PolFileToMPRealVector(mps_context *s, mps_monomial_poly *p, slong wp) {
  mpfr::mpreal::set_default_prec(3000); // TODO(orebas): magic number
  rdpe_t u;
  // cdpe_t cx;
  std::vector<std::complex<mpfr::mpreal>> results;
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
    results.push_back(coeff);
  }
  mpfr_clear(mpfrre);
  mpfr_clear(mpfrim);
}