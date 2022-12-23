#include <cassert>
#include <random>

#include "acb_calc.h"
#include "acb_poly.h"
#include "arb.h"
#include "arbxx.hpp"
#include <iostream>
#include <iterator>
#include <map>

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
                  calcCoeff(p - 1,
                            2 * k - j); // TODO use better function to times 2
            } else {
              result -= calcCoeff(p - 1, j) * calcCoeff(p - 1, 2 * k - j);
            }
          }
        } else {
          for (slong j = k + 1; j <= 2 * k; j++) {
            if (j % 2 == 0) {
              result +=
                  calcCoeff(p - 1, j) *
                  calcCoeff(p - 1,
                            2 * k - j); // TODO use better function to times 2
            } else {
              result -= calcCoeff(p - 1, j) * calcCoeff(p - 1, 2 * k - j);
            }
          }
        }
        acb_mul_2exp_si(result.c, result.c, 1); // TODO make a member function.
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
             i++) { // TODO this is too many iteration, trim it down
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

      acb_poly_shift_right(t3.pol, t3.pol, 1); // TODO: "unsquare x^2"

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
        outp[i].push_back(absval); // TODO FIX ARB
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

template <typename T>
std::ostream &operator<<(std::ostream &out, const std::vector<T> &v) {
  if (!v.empty()) {
    out << '[';
    std::copy(v.begin(), v.end(), std::ostream_iterator<T>(out, ", "));
    out << "\b\b]";
  }
  out << std::endl;
  return out;
}

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

std::vector<ARB> ACBVectorComp(
    std::vector<ACB> &a,
    std::vector<ACB>
        &b) { // compare b with a as the reference  TODO make this const later
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

  for (int i = 0; i < a.size(); i++) {
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
    std::random_device
        rd; // TODO make this a static object and only have one generator
    std::mt19937 mt(rd());
    std::uniform_real_distribution<double> dist(-1.0, 1.0);
    std::vector<ACB> roots;
    // TODO test strip initial zeros
    // TODO make monic
    // TODO replace x^n by x if you can (i.e. even polynomials)

    slong stripped_zeros = acb_poly_valuation(poly.pol);
    if (stripped_zeros > 0) {
      for (int i = 0; i < stripped_zeros; i++) {
        roots.push_back(ACB(0.0, 0.0, poly.intprec));
      }
      acb_poly_shift_right(poly.pol, poly.pol, stripped_zeros);
    }
    auto maxRad = RRE.maxRadius(depth); // TODO 20 is a magic number
    auto xrand = dist(rd);
    auto yrand = dist(rd);
    std::cout << xrand << yrand;
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
