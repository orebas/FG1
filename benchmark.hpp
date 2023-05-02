
#pragma once
#include <nlohmann/json.hpp>
using json = nlohmann::json;

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

#include <cassert>
#include <random>

#include "acb_calc.h"
#include "acb_poly.h"
#include "arb.h"
#include "arbxx.hpp"
#include <iostream>
#include <iterator>
#include <map>

template <class T>
std::ostream &operator<<(std::ostream &out, const std::complex<T> &v) {
  out << "[" << v.real() << ", " << v.imag() << "]" << std::endl;
  return out;
}

template <class T, class V>
std::ostream &operator<<(std::ostream &out, const std::pair<T, V> &v) {
  out << "[" << v.first << ", " << v.second << "]" << std::endl;
  return out;
}

template <typename T>
std::ostream &operator<<(std::ostream &out, const std::vector<T> &v) {
  if (!v.empty()) {
    out << '[';
    for (int i = 0; i < v.size() - 1; i++) {
      out << v[i] << ", ";
    }
    out << v[v.size() - 1];
    out << "]";
  }
  out << std::endl;
  return out;
}

template <typename TimeT = std::chrono::milliseconds> struct measure {
  template <typename F, typename... Args>
  static typename TimeT::rep execution(F func, Args &&...args) {
    auto start = std::chrono::steady_clock::now();

    // Now call the function with all the parameters you need.
    func(std::forward<Args>(args)...);

    auto duration = std::chrono::duration_cast<TimeT>(
        std::chrono::steady_clock::now() - start);

    return duration.count();
  }

  template <typename F, typename... Args>
  static typename TimeT::rep execution2(F func, Args &&...args) {
    auto start = std::chrono::steady_clock::now();

    // Now call the function with all the parameters you need.
    func(std::forward<Args>(args)...);

    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(
        std::chrono::steady_clock::now() - start);

    auto start2 = std::chrono::steady_clock::now();
    long iters = 1000.0 / (duration.count() + 1);
    iters = std::max(1L, iters);
    iters = std::min(iters, 50L);
    for (long i = 0; i < iters; ++i) {
      func(std::forward<Args>(args)...);
    }
    auto duration2 = std::chrono::duration_cast<TimeT>(
        std::chrono::steady_clock::now() - start2);
    return duration2.count() / iters;
  }
};

void usage(mps_context *s, const char *program);
void parsePol(const std::string &polfilename);
typedef mps_context *mpscp;

void *cleanup_context(mpscp ctx, void *user_data, mps_polynomial *poly,
                      mpscp &s, bool outp);

void DLGFG(const ComplexPoly &polarg, std::vector<ComplexPoly> pvec,
           std::vector<ComplexPoly> qvec, long depth);
void solveDebug(mps_context *local_s, mps_polynomial *local_poly,
                const std::string &polfilename);
void solveCompare(mps_context *local_s, mps_polynomial *local_poly,
                  const std::string &polfilename);

class RecursiveRadiiEstimator {
public:
  ComplexPoly p0;
  ComplexPoly q0;
  slong prec;
  slong degree;
  std::map<std::pair<slong, slong>, ACB> coefficients,
      qcoefficients; // stores DLG coefficients, if and only if they have been
                     // calculated.
  // calcCoff calculates the coefficient of x^i in p_k (where p_0 is initial
  // polynomial), stores it in coefficients, and returns it.

  ACB calcCoeff(slong p, slong k) {
    ACB result(0, 0, prec);
    std::pair<slong, slong> key(p, k);
    if (k >= 0 && k <= degree) {
      if (p == 0) {
        result = p0.getCoeff(k);
      } else if (coefficients.contains(key)) {
        result = coefficients.at(key); // because ACB not default constructible.
      } else {
        if (2 * k > degree) {
          for (slong j = 0; j <= k - 1; j++) {
            if (j % 2 == 0) {
              result +=
                  calcCoeff(p - 1, j) *
                  calcCoeff(
                      p - 1,
                      2 * k - j); // TODO(orebas) use better function to times 2
            } else {
              result -= calcCoeff(p - 1, j) * calcCoeff(p - 1, 2 * k - j);
            }
          }
        } else {
          for (slong j = k + 1; j <= 2 * k; j++) {
            if (j % 2 == 0) {
              result +=
                  calcCoeff(p - 1, j) *
                  calcCoeff(
                      p - 1,
                      2 * k - j); // TODO(orebas) use better function to times 2
            } else {
              result -= calcCoeff(p - 1, j) * calcCoeff(p - 1, 2 * k - j);
            }
          }
        }
        acb_mul_2exp_si(result.c, result.c,
                        1); // TODO(orebas) make a member function.
        if (k % 2 == 0) {
          result += calcCoeff(p - 1, k).square();
        } else {
          result -= calcCoeff(p - 1, k).square();
        }
        coefficients.insert_or_assign(key, result);
      }
    }

    return result;
  }

  ACB calcQcoeff(slong q, slong k) {
    ACB result(0, 0, prec);
    std::pair<slong, slong> key(q, k);
    if (k >= 0 && k <= degree) {
      if (q == 0) {
        result = q0.getCoeff(k);
      } else if (qcoefficients.contains(key)) {
        result = qcoefficients.at(key);
      } else {
        for (slong i = 0; i <= 2 * k + 1;
             i++) { // TODO(orebas) this is too many iteration, trim it down
          slong j = 2 * k + 1 - i;

          if (i % 2 == 0) {
            result += calcQcoeff(q - 1, i) * calcCoeff(q - 1, j);
          } else {
            result -= calcQcoeff(q - 1, i) * calcCoeff(q - 1, j);
          }
          if (j % 2 == 0) {
            result -= calcQcoeff(q - 1, i) * calcCoeff(q - 1, j);
          } else {
            result += calcQcoeff(q - 1, i) * calcCoeff(q - 1, j);
          }
        }
        acb_mul_2exp_si(result.c, result.c, -1);
        qcoefficients.insert_or_assign(key, result);
      }
    }

    return result;
  }

  ARB minRadius(slong depth) {
    ARB r(0, prec);
    r = (calcCoeff(depth, 0) / calcCoeff(depth, 1))
            .abs()
            .root_ui(std::pow(2, depth));
    return r;
  }

  ARB maxRadius(slong depth) {
    ARB r(0, prec);
    r = (calcCoeff(depth, degree - 1) / calcCoeff(depth, degree))
            .abs()
            .root_ui(std::pow(2, depth));
    return r;
  }

  ACB minRoot(slong depth) {
    ACB y0 = calcQcoeff(depth, 0);
    ACB y1 = calcCoeff(depth, 1);
    ACB root_est = y0 / y1;
    return root_est;
  }

  ACB maxRoot(slong depth) {
    ACB y0 = calcQcoeff(depth, degree - 1);
    ACB y1 = calcCoeff(depth, degree);
    ACB z0 = calcQcoeff(depth, degree - 2);
    ACB z1 = calcCoeff(depth, degree - 1);
    ACB t = y0 / y1 - z0 / z1;
    return t;
  }

  RecursiveRadiiEstimator(const ComplexPoly &polargtocopy)
      : p0(polargtocopy), q0(polargtocopy), prec(polargtocopy.intprec),
        degree(p0.degree()) {
    acb_poly_derivative(q0.pol, q0.pol, prec);
    acb_poly_shift_left(q0.pol, q0.pol, 1);
  }
};

class FGSolver {
public:
  ComplexPoly p0;
  std::vector<ComplexPoly> pvec;
  std::vector<ComplexPoly> qvec;
  slong depth, prec;

public:
  void generatePnQnPolynomials() { // depth and prec must already be set
    ComplexPoly negx(prec);
    acb_poly_set_coeff_si(negx.pol, 1, -1);
    for (slong i = 1; i < depth; i++) {
      // std::cout << "\n\ni is " << i << std::endl;
      pvec[i].explicitCopy(pvec[i - 1].graeffe());
      // pvec[i].print();
    }
    for (slong i = 1; i < depth; i++) {
      ComplexPoly t1(prec), t2(prec), t3(prec);
      ComplexPoly qnegx = qvec[i - 1].compose(negx);
      ComplexPoly pnegx = pvec[i - 1].compose(negx);

      // acb_poly_compose(qnegx.pol, qvec[i - 1].pol, negx.pol, qvec[i -
      // 1].intprec); acb_poly_compose(pnegx.pol, pvec[i - 1].pol, negx.pol,
      // pvec[i - 1].intprec);
      t1 = qvec[i - 1] * pnegx;
      // acb_poly_mul(t1.pol, qvec[i - 1].pol, pnegx.pol, qvec[i - 1].intprec);
      // acb_poly_mul(t2.pol, qnegx.pol, pvec[i - 1].pol, pvec[i - 1].intprec);
      t2 = qnegx * pvec[i - 1];

      // acb_poly_sub(t3.pol, t1.pol, t2.pol, t1.intprec);
      t3 = t1 - t2;
      // acb_poly_neg(t3.pol, t3.pol);
      // std::cout << "T3A     " << t3.intprec << "     ";
      // t3.print(3);

      acb_poly_scalar_mul_2exp_si(
          t3.pol, t3.pol, -1); // at this point t3 should be divisible by x

      acb_poly_shift_right(t3.pol, t3.pol, 1); // TODO(orebas): "unsquare x^2"

      t3.unsquare();

      qvec[i].explicitCopy(t3);
    }
  }
  FGSolver(const ComplexPoly &polargtocopy, long initdepth)
      : p0(polargtocopy), depth(initdepth), prec(polargtocopy.intprec) {
    pvec.resize(depth, ComplexPoly(prec));
    qvec.resize(depth, ComplexPoly(prec));
    pvec[0].explicitCopy(polargtocopy);
    qvec[0].explicitCopy(polargtocopy);
    acb_poly_derivative(qvec[0].pol, qvec[0].pol, qvec[0].intprec);
    acb_poly_shift_left(qvec[0].pol, qvec[0].pol, 1); // q_0 = xp'(x)
  }
  void DLGRadii(std::vector<std::vector<ARB>> &outp) {
    for (slong i = 1; i < depth; i++) {
      for (int j = 0; j <= acb_poly_degree(pvec[i].pol) - 1; j++) {
        ACB y0(0, 0, prec), y1(0, 0, prec);
        ARB absval(0, prec);

        y0 = pvec[i].getCoeff(j);
        y1 = pvec[i].getCoeff(j + 1);

        acb_div(y1.c, y0.c, y1.c, y0.intprec);
        absval = y1.abs();
        arb_root_ui(absval.r, absval.r, std::pow(2, i), y1.intprec);
        outp[i].push_back(absval); // TODO(orebas) FIX ARB
                                   // std::cout << "[" << j << "]  ";
                                   // arb_printd(absval, 30);
                                   // std::cout << " " << std::endl;
      }
    }
    // std::cout << std::endl;
  }
  void FGRoots(std::vector<std::vector<ACB>> &outp) {
    for (slong k = 1; k < depth; k++) {
      ACB y0(0, 0, prec), y1(0, 0, prec), prev_rat(0, 0, prec);
      ACB root_est(0, 0, prec);
      acb_poly_get_coeff_acb(y0.c, qvec[k].pol, 0);
      acb_poly_get_coeff_acb(y1.c, pvec[k].pol, 1);
      acb_div(root_est.c, y0.c, y1.c, root_est.intprec);
      acb_set(prev_rat.c, root_est.c);
      // std::cout << "[" << 0 << "]  ";
      // root_est.print(20);
      outp[k].push_back(root_est);
      // std::cout << " " << std::endl;

      for (int j = 1; j < acb_poly_degree(qvec[0].pol); j++) {
        acb_poly_get_coeff_acb(y0.c, qvec[k].pol, j);
        acb_poly_get_coeff_acb(y1.c, pvec[k].pol, j + 1);
        acb_div(y0.c, y0.c, y1.c, y0.intprec);
        acb_sub(root_est.c, y0.c, prev_rat.c, root_est.intprec);
        acb_set(prev_rat.c, y0.c);
        // std::cout << "[" << j << "]  ";
        // root_est.print(20);
        // std::cout << " " << std::endl;
        outp[k].push_back(root_est);
      }
    }
  }
};

/*template <typename T>
std::ostream &operator<<(std::ostream &out, const std::vector<T> &v) {
  if (!v.empty()) {
    out << '[';
    std::copy(v.begin(), v.end(), std::ostream_iterator<T>(out, ", "));
    out << "\b\b]";
  }
  out << std::endl;
  return out;
}*/  //moved to jsonpoly.hpp

template <typename T, size_t s>
std::ostream &operator<<(std::ostream &out, const std::array<T, s> &v) {
  if (!v.empty()) {
    out << '[';
    std::copy(v.begin(), v.end(), std::ostream_iterator<T>(out, ", "));
    out << "\b\b]";
  }
  out << std::endl;
  return out;
}

std::ostream &operator<<(std::ostream &ost, const ARB &ls) {
  auto a = arb_get_str(ls.r, 20, 0);
  ost << a << " ";
  free(a);
  return ost;
}

std::ostream &operator<<(std::ostream &ost, const ACB &ls) {
  auto a = arb_get_str(acb_realref(ls.c), 20, 0);
  auto b = arb_get_str(acb_imagref(ls.c), 20, 0);
  ost << a << " " << b << " " << std::endl;
  free(a);
  free(b);
  return ost;
}

std::vector<ARB>
ACBVectorComp(std::vector<ACB> &a,
              std::vector<ACB> &b) { // compare b with a as the reference
                                     // TODO(orebas) make this const later
  std::vector<ARB> results;
  slong prec = a[0].intprec;
  // std::cout << a.size() << " " << b.size() << std::endl;
  if (a.size() != b.size()) {
    std::cout << "Different Sizes!!!" << std::endl;
    return results;
  }
  std::vector<ARB> absdiff(a.size(), ARB(0, prec));
  std::vector<ARB> rel(a.size(), ARB(0, prec));

  ARB maxrel(0, prec), maxabs(0, prec), rmsrel(0, prec), rmsabs(0, prec);

  for (std::size_t i = 0; i < a.size(); i++) {
    // std::cout << "A:" << a[i] << "b:" << b[i] << " ";
    auto temp1 = a[i] - b[i];
    auto temp2 = temp1.abs();
    // std::cout << "absdiff temp1" << temp1 << "absdiff 2:" << temp2;
    absdiff[i] = (a[i] - b[i]).abs();
    rel[i] = absdiff[i] / (a[i].abs());
    // std::cout << "abs:" << absdiff[i] << " res:" << rel[i];
    if (arb_ge(rel[i].r, maxrel.r)) {
      maxrel = rel[i];
    }
    if (arb_ge(absdiff[i].r, maxabs.r)) {
      maxabs = absdiff[i];
    }
    rmsrel += rel[i] * rel[i];
    rmsabs += absdiff[i] * absdiff[i];
  }
  rmsrel /= ARB((double)a.size(), prec);
  rmsabs /= ARB((double)a.size(), prec);
  results.push_back(maxabs);
  results.push_back(maxrel);
  results.push_back(rmsabs);
  results.push_back(rmsrel);

  rmsrel = rmsrel.root_ui(2);
  rmsabs = rmsabs.root_ui(2);

  return results;
}

class RootSolver1 {
public:
  ComplexPoly poly;
  slong depth;

  RecursiveRadiiEstimator RRE;
  FGSolver FGS;
  std::vector<std::vector<ACB>> rootapprox;
  RootSolver1(const ComplexPoly &intpoly, slong initdepth)
      : poly(intpoly), depth(initdepth), RRE(poly), FGS(poly, initdepth),
        rootapprox(depth) {}

  std::vector<ACB> solve(void) {
    std::random_device rd; // TODO(orebas) make this a static object and only
                           // have one generator
    std::mt19937 mt(rd());
    std::uniform_real_distribution<double> dist(-1.0, 1.0);
    std::vector<ACB> roots;
    // TODO(orebas) test strip initial zeros
    // TODO(orebas) make monic
    // TODO(orebas) replace x^n by x if you can (i.e. even polynomials)

    slong stripped_zeros = acb_poly_valuation(poly.pol);
    if (stripped_zeros > 0) {
      for (int i = 0; i < stripped_zeros; i++) {
        roots.push_back(ACB(0.0, 0.0, poly.intprec));
      }
      acb_poly_shift_right(poly.pol, poly.pol, stripped_zeros);
    }
    auto maxRad = RRE.maxRadius(depth); // TODO(orebas) 20 is a magic number
    auto xrand = dist(rd);
    auto yrand = dist(rd);
    //std::cout << xrand << yrand;
    ACB randomcenter = ACB(maxRad) * ACB(xrand, yrand, poly.intprec);
    ComplexPoly shiftedpoly(poly);
    acb_poly_taylor_shift(shiftedpoly.pol, shiftedpoly.pol, randomcenter.c,
                          shiftedpoly.intprec);
    FGSolver FGS(shiftedpoly, depth);
    FGS.generatePnQnPolynomials();
    FGS.FGRoots(rootapprox);
    for (auto &v : rootapprox) {
      for (auto &i : v) {
        i += randomcenter;
      }
    }
    for (auto i : rootapprox.back()) {
      roots.push_back(i);
    }

    return roots;
  }
};

class SimpleNewtonSolver {
  std::random_device rd; //
  std::mt19937 mt;
  std::uniform_real_distribution<double> dist;

public:
  SimpleNewtonSolver() : mt(rd()), dist(-1.0, 1.0) {}

  std::vector<ACB> solve(ComplexPoly &poly, slong prec) {
    assert(prec > 20);
    const slong n = poly.degree();
    RecursiveRadiiEstimator RRE(poly);
    auto maxRad = ACB(RRE.maxRadius(10)); // TODO(orebas) MAGIC NUMBER
    // std::cout << "Radius estimate is: ";
    // maxRad.print();
    // std::cout << std::endl;

    // mag_t fujibound;
    // mag_init(fujibound);
    // acb_poly_root_bound_fujiwara(fujibound, poly.pol);
    // ARB fbound(0, prec);
    // arb_set_interval_mag(fbound.r, fujibound, fujibound, prec);
    // auto xrand = dist(rd);
    // auto yrand = dist(rd);
    // ACB randomcenter = ACB(fbound) * ACB(xrand, yrand, poly.intprec);
    std::vector<ACB> roots;
    ACB pval = (ACB(0, 0, prec));
    ACB pprimeval = pval;
    ACB newpval = pval;
    ACB x = pval;
    ACB tempsum = x;
    ACB one(1, 0, prec);
    ACB eps(1, 0, prec);
    ACB sqeps(1, 0, prec);
    ACB dir(0, 0, prec);
    ACB dir2(0, 0, prec);
    acb_mul_2exp_si(eps.c, one.c, -prec);
    acb_mul_2exp_si(sqeps.c, one.c, -((prec / 2) + 1)); // integer division
    ACB twoeps(eps);
    twoeps *= ACB(2, 0, prec);
    acb_one(one.c);
    for (slong k = 0; k < n; k++) {
      x = maxRad * ACB(dist(rd), dist(rd), poly.intprec); // random start!

      for (slong i = 0; i < (n * 2) + 100; i++) {
        // we do maximum n*2+50 iterations, but
        // break earlier if we need to.
        acb_poly_evaluate2(pval.c, pprimeval.c, poly.pol, x.c, prec);
        if (pval.abs() < (eps * twoeps).abs()) {
          break;
        }
        if (pprimeval.abs() < sqeps.abs()) {
          acb_sgn(pprimeval.c, pprimeval.c, prec);
        }
        if (pprimeval.abs() < eps.abs()) {
          pprimeval = sqeps * ACB(dist(rd), dist(rd), poly.intprec); // random?
        }
        dir = pval / pprimeval;

        acb_zero(tempsum.c);
        for (std::size_t j = 0; j < roots.size(); j++) {
          if ((x - roots[j]).abs() > sqeps.abs()) {
            tempsum += one / (x - roots[j]);
          }
        }
        dir2 = dir * (one / (one - (dir * tempsum)));
        x -= dir2;
        acb_get_mid(x.c, x.c);
      }
      acb_poly_evaluate(newpval.c, poly.pol, x.c, prec);
      if (newpval.abs() < pval.abs()) {
        roots.push_back(x);
      } else {
        roots.push_back(x + dir2);
      }
    }
    //  mag_clear(fujibound);
    sortRootVector(roots);
    // for (slong k = 0; k < roots.size(); k++) {
    //   roots[k].print(20);
    //   std::cout << std::endl;
    // }
    // std::cout << std::endl;

    return roots;
  }
};

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
  //std::cout << radii;
  //std::cout << rootapprox;
  // std::cout << powersums(110, 20, 6000);
  //std::cout << "Begin RRE debugging" << std::endl;
  RecursiveRadiiEstimator RRE(acb_style_poly);

  /*for (int d = 0; d < 15; d++) {
      std::cout << "RRE " << d << " " << std::endl;
      FGS.pvec[d].print();
      for (int n = 0; n <= acb_style_poly.degree(); n++) {
          std::cout << RRE.calcCoeff(d, n);
      }
  }*/

  for (int depth = 0; depth < 15; depth++) {
    //std::cout << depth << "[min]" << RRE.minRadius(depth);
    //std::cout << "[max]" << RRE.maxRadius(depth);
  }

  for (int depth = 0; depth < 15; depth++) {
    //std::cout << depth << "[minRoot]" << RRE.minRoot(depth);

    //std::cout << depth << "[minRoot]" << RRE.maxRoot(depth);
  }

  double initTime = measure<std::chrono::milliseconds>::execution(func);
  (void)fflush(stdout);
  //std::cout << polfilename << " took " << initTime << "ms." << std::endl;
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

void timeSolvers(std::vector<ComplexPoly> pvec);

inline std::vector<ACB> powersums(slong depth, slong max_k, slong prec) {
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