#pragma once

#include "acb_calc.h"
#include "acb_poly.h"
#include "arb.h"
#include <cassert>
template <class T = double> class Array3d {
private:
  std::size_t w;
  std::size_t h;
  std::size_t d;

public:
  std::vector<T> data;

  Array3d(std::size_t width, std::size_t height, std::size_t depth)
      : w(width), h(height), d(depth), data(w * h * d, 0) {}

  Array3d(std::size_t width, std::size_t height, std::size_t depth,
          const T &initializer)
      : w(width), h(height), d(depth), data(w * h * d, initializer) {}

  inline T &at(std::size_t x, std::size_t y, std::size_t z) {
    return data[x * h * d + y * d + z];
  }

  inline T at(std::size_t x, std::size_t y, std::size_t z) const {
    return data[x * h * d + y * d + z];
  }

  inline T &operator()(std::size_t x, std::size_t y, std::size_t z) {
    return data[x * h * d + y * d + z];
  }

  inline T operator()(std::size_t x, std::size_t y, std::size_t z) const {
    return data[x * h * d + y * d + z];
  }
  inline std::size_t width() { return w; }
  inline std::size_t height() { return h; }
  inline std::size_t depth() { return d; }
};

template <class T = double> class Array4d {
private:
  std::size_t w;
  std::size_t h;
  std::size_t d;
  std::size_t l;

public:
  std::vector<T> data;

  Array4d(std::size_t r, std::size_t s, std::size_t t, std::size_t u)
      : w(r), h(s), d(t), l(u), data(w * h * d * l, 0) {}

  inline T &at(std::size_t x, std::size_t y, std::size_t z, std::size_t t) {
    return data[x * h * d * l + y * d * l + z * l + t];
  }

  inline T at(std::size_t x, std::size_t y, std::size_t z,
              std::size_t t) const {
    return data[x * h * d * l + y * d * l + z * l + t];
  }

  inline T &operator()(std::size_t x, std::size_t y, std::size_t z,
                       std::size_t t) {
    return data[x * h * d * l + y * d * l + z * l + t];
  }

  inline T operator()(std::size_t x, std::size_t y, std::size_t z,
                      std::size_t t) const {
    return data[x * h * d * l + y * d * l + z * l + t];
  }

  inline std::size_t internalRef(std::size_t x, std::size_t y, std::size_t z,
                                 std::size_t t) {
    return x * h * d * l + y * d * l + z * l + t;
  }
};

class ARF {
public:
  arf_t f;
  slong intprec = 1;
  inline ::arf_ptr arf_ptr_int() { return f; }
  ARF(const arf_t &s, slong ip) : intprec(ip) {
    arf_init(f);
    arf_set(f, s);
  }

  ARF(double x, slong ip) : intprec(ip) {
    arf_init(f);
    arf_set_d(f, x);
  }

  ARF(const ARF &other) : ARF(other.f, other.intprec) {}
  ARF(ARF &&other) noexcept : intprec(other.intprec) {
    arf_swap(arf_ptr_int(), other.arf_ptr_int());
    other.intprec = -1;
  }

  ARF &operator=(const ARF &other) { return *this = ARF(other); }

  ARF &operator=(ARF &&other) noexcept {
    // std::swap(c, other.c);
    arf_swap(arf_ptr_int(),
             other.arf_ptr_int()); // Do I need to deinitialize here?
    return *this;
  }

  ~ARF() {
    if (intprec > 0) {
      arf_clear(f);
    }
  }
};

class ARB {
public:
  arb_t r;
  slong intprec = 1;

  inline ::arb_ptr arb_ptr_int() { return r; }

  ARB(const arb_t &s, slong ip) : intprec(ip) {
    arb_init(r);
    arb_set(r, s);
  }
  ARB(double x, slong ip) : intprec(ip) {
    arb_init(r);
    arb_set_d(r, x);
  }
  ARB(int x, slong ip) : intprec(ip) {
    arb_init(r);
    arb_set_d(r, x);
  }
  ARB(const ARB &other) : ARB(other.r, other.intprec) {}
  ARB(ARB &&other) noexcept : intprec(other.intprec) {
    arb_swap(arb_ptr_int(), other.arb_ptr_int());
    other.intprec = -1;
  }
  ARB(const arf_t arf, slong prec) : intprec(prec) {
    arb_init(r);
    arb_set_arf(r, arf);
  }
  ARB(const mpfr_t mpfr, slong prec) : intprec(prec) {
    arb_init(r);
    arb_set_interval_mpfr(r, mpfr, mpfr, prec);
  }

  ARB &operator=(const ARB &other) { return *this = ARB(other); }

  ARB root_ui(slong e) {
    ARB temp(0.0, intprec);
    arb_root_ui(temp.r, this->r, e, intprec);
    return temp;
  }

  ARB &operator=(ARB &&other) noexcept {
    // std::swap(c, other.c);
    intprec = other.intprec;
    arb_swap(arb_ptr_int(),
             other.arb_ptr_int()); // Do I need to deinitialize here?
    other.intprec = -1;
    return *this;
  }

  ARB operator/=(const ARB &other) {
    arb_div(this->r, this->r, other.r, intprec);
    return *this;
  }

  ARB operator/(const ARB &other) { return ARB(*this) /= other; }
  ARB operator-=(const ARB &other) {
    arb_sub(this->r, this->r, other.r, intprec);
    return *this;
  }

  ARB operator-(const ARB &other) { return ARB(*this) -= other; }

  ARB operator+=(const ARB &other) {
    arb_add(this->r, this->r, other.r, intprec);
    return *this;
  }

  ARB operator+(const ARB &other) { return ARB(*this) += other; }

  ARB operator*=(const ARB &other) {
    arb_mul(this->r, this->r, other.r, intprec);
    return *this;
  }
  ARB operator*(const ARB &other) { return ARB(*this) *= other; }

  bool operator<(const ARB &other) { return (arb_lt(this->r, other.r) != 0); }
  bool operator>(const ARB &other) { return (arb_gt(this->r, other.r) != 0); }

  ~ARB() {
    if (intprec > 0) {
      arb_clear(r);
    }
  }
  void explicitCopy(const ARB &a) {
    arb_set(r, a.r);
    intprec = a.intprec;
  }
  void print() const { arb_printd(r, 30); }
};

ARB inline fabs(ARB &a) {
  ARB result(a);
  ARB b(a);
  arb_neg(b.r, b.r);
  if (a < b) {
    result = b;
  } else {
    result = a;
  }
  return result;
}

class ACB {
public:
  acb_t c;
  slong intprec = 1;

  inline ::acb_ptr acb_ptr_int() { return c; }
  // inline ::mpfr_srcptr mpreal::mpfr_ptr() const { return mp; }
  // inline ::mpfr_srcptr mpreal::mpfr_srcptr() const { return mp; }

  /*ACB(slong ip) : intprec(ip) {
      acb_init(c);  // this already sets to zero.
  }*/
  ACB(const acb_t &s, slong ip) : intprec(ip) {
    acb_init(c);
    acb_set(c, s);
  }
  ACB(double x, double y, slong ip) : intprec(ip) {
    acb_init(c);
    acb_set_d_d(c, x, y);
  }

  ACB(const ARB &r) : intprec(r.intprec) {
    acb_init(c);
    acb_set_arb(this->c, r.r);
  }

  ACB(const ARB &r, const ARB &i) : intprec(r.intprec) {
    acb_init(c);
    acb_set_arb_arb(this->c, r.r, i.r);
  }

  /* ACB(const mpc_t &cold, slong prec) : intprec(prec)
  {
      mpfr_t rempfr, immpfr;
      arf_t rearf, imarf;
      // arf_init(rearf);
      // arf_init(imarf);
      mpfr_init(rempfr);
      mpfr_init(immpfr);
      mpfr_set_f(rempfr, cold->r);
      mpfr_set_f(immpfr, cold->i);
      // arf_set_mpfr(rearf, rempfr);
      // arf_set_mpfr(imarf, immpfr);

      ACB result(ARB(rempfr, prec), (ARB(immpfr, prec)));
      mpfr_clear(rempfr);
      mpfr_clear(immpfr);
      // arf_clear(rearf);
      // arf_clear(imarf);
      return result;
  }
  */

  ACB(int a, int b, slong ip) : intprec(ip) {
    acb_init(c);
    acb_set_si_si(c, a, b);
  }

  /* ACB(slong a, slong b, slong ip) : intprec(ip)
   {
       acb_init(c);
       acb_set_si_si(c, a, b);
   }*/
  ACB(const mpc_t &p, slong local_prec) : intprec(local_prec) {
    acb_init(c);
    mpfr_t mpfrre;
    mpfr_t mpfrim;
    arf_t arfre;
    arf_t arfim;
    mpfr_init2(mpfrre, local_prec);
    mpfr_init2(mpfrim, local_prec);
    arf_init(arfre);
    arf_init(arfim);
    ARB arbre(0, local_prec);
    ARB arbim(0, local_prec);
    mpfr_set_f(mpfrre, (p)->r, MPFR_RNDN);
    mpfr_set_f(mpfrim, (p)->i, MPFR_RNDN);
    arf_set_mpfr(arfre, mpfrre);
    arf_set_mpfr(arfim, mpfrim);
    arb_set_arf(arbre.r, arfre);
    arb_set_arf(arbim.r, arfim);
    acb_set_arb_arb(c, arbre.r, arbim.r);
    mpfr_clear(mpfrre);
    mpfr_clear(mpfrim);
    arf_clear(arfre);
    arf_clear(arfim);
  }

  ACB(const ACB &other) : ACB(other.c, other.intprec) {}
  ACB(ACB &&other) noexcept : intprec(other.intprec) {
    // acb_ptr_int() = 0;
    acb_swap(acb_ptr_int(), other.acb_ptr_int());
    other.intprec = -1;
    // c = other.c;
    // other.c = nullptr;
  }

  ACB &operator=(const ACB &other) { return *this = ACB(other); }

  ACB &operator=(ACB &&other) noexcept {
    // std::swap(c, other.c);
    intprec = other.intprec;
    acb_swap(acb_ptr_int(), other.acb_ptr_int());
    other.intprec = -1;
    return *this;
  }
  ACB operator+=(const ACB &other) {
    acb_add(this->c, this->c, other.c, intprec);
    return *this;
  }
  ACB operator-=(const ACB &other) {
    acb_sub(this->c, this->c, other.c, intprec);
    return *this;
  }
  ACB operator/=(const ACB &other) {
    acb_div(this->c, this->c, other.c, intprec);
    return *this;
  }
  ACB square() {
    acb_sqr(this->c, this->c, intprec);
    return *this;
  }
  ACB operator*=(const ACB &other) {
    acb_mul(this->c, this->c, other.c, intprec);
    return *this;
  }
  ACB operator*(const ACB &other) { return ACB(*this) *= other; }
  ACB operator+(const ACB &other) { return ACB(*this) += other; }
  ACB operator/(const ACB &other) { return ACB(*this) /= other; }
  ACB operator-(const ACB &other) { return ACB(*this) -= other; }
  ACB power(slong p) {
    ACB temp(0, 0, intprec);
    acb_pow_si(temp.c, this->c, p, intprec);
    return temp;
  }
  ARB abs() const {
    assert(this->intprec > 0);
    ARB temp(0, intprec);
    // std::cout << "temp init:" << arb_get_str(temp.r, 20, 0) << std::endl;
    // std::cout << "cplx init:" << arb_get_str(acb_realref(this->c), 20, 0) <<
    // " " << arb_get_str(acb_imagref(this->c), 20, 0);
    acb_abs(temp.r, this->c, this->intprec);
    // std::cout << "precision:" << this->intprec << std::endl;
    // std::cout << std::endl
    //<< "abs result:" << arb_get_str(temp.r, 20, 0) << std::endl;

    return temp;
  }

  ARF abs_ubound_arf() const {
    ARF l(0, intprec);
    acb_get_abs_ubound_arf(l.f, this->c, intprec);
    return l;
  }

  ~ACB() {
    if (intprec > 0) {
      acb_clear(c);
    }
  }

  void explicitCopy(const ACB &a) { acb_set(c, a.c); }
  void print(int digits = 5) { acb_printd(c, digits); }
};

// ACB operator*(const ACB &lhs, const ACB &rhs) {
//     return ACB(lhs) *= rhs;
// }

void inline sortRootVector(std::vector<ACB> &roots) {
  struct {
    bool operator()(const ACB &a, const ACB &b) const {

      ARB re1(0.0, a.intprec);
      ARB im1(0.0, b.intprec);

      ARB re2(0.0, a.intprec);
      ARB im2(0.0, b.intprec);

      acb_get_real(re1.r, a.c);
      acb_get_real(re2.r, b.c);
      acb_get_imag(im1.r, a.c);
      acb_get_imag(im2.r, b.c);

      if (arb_lt(re1.r, re2.r) != 0) {
        return true;
      }
      if (arb_gt(re1.r, re2.r) != 0) {
        return false;
      }
      if (arb_lt(im1.r, im2.r) != 0) {
        return true;
      }
      if (arb_gt(im1.r, im2.r) != 0) {
        return false;
      }
      return false;
    }
  } customLessLex;

  std::sort(roots.begin(), roots.end(), customLessLex);
}

class ComplexPoly {
public:
  acb_poly_t pol;
  slong intprec = 1;

  typedef acb_poly_struct *acb_poly_ptr;

  acb_poly_ptr acb_poly_ptr_int() { return pol; }

  explicit ComplexPoly(slong ip) : intprec(ip) { acb_poly_init(pol); }
  ComplexPoly(const ComplexPoly &other) : intprec(other.intprec) {
    acb_poly_init(pol);
    acb_poly_set(pol, other.pol);
  }
  ComplexPoly(ComplexPoly &&other) noexcept : intprec(other.intprec) {
    acb_poly_swap(acb_poly_ptr_int(), other.acb_poly_ptr_int());
    other.intprec = -1;
  }
  ComplexPoly &operator=(const ComplexPoly &other) {
    return *this = ComplexPoly(other);
  }
  ComplexPoly &operator=(ComplexPoly &&other) noexcept {
    intprec = other.intprec;
    acb_poly_swap(acb_poly_ptr_int(), other.acb_poly_ptr_int());
    other.intprec = -1;
    return *this;
  }

  ComplexPoly compose(const ComplexPoly &polarg) {
    ComplexPoly result(intprec);
    acb_poly_compose(result.pol, this->pol, polarg.pol, intprec);
    return result;
  }
  ComplexPoly(mps_context *s, mps_monomial_poly *p, std::size_t wp)
      : intprec(wp) {
    acb_poly_init(pol);
    // int i;
    // rdpe_t apol, ax
    rdpe_t u;
    // cdpe_t cx;

    pthread_mutex_lock(&p->mfpc_mutex[0]);
    if (mpc_get_prec(p->mfpc[0]) < wp) {
      // DEBUG
      // std::cout << "debug :" << p->density << std::endl;
      // END DEBUG
      std::cout << "it is " << mpc_get_prec(p->mfpc[0]) << " and " << wp
                << std::endl;
      pthread_mutex_unlock(&p->mfpc_mutex[0]);
      mps_monomial_poly_raise_precision(s, MPS_POLYNOMIAL(p), wp);
      std::cout << "it is " << mpc_get_prec(p->mfpc[0]) << " and " << wp
                << std::endl;
    } else
      pthread_mutex_unlock(&p->mfpc_mutex[0]);

    // if (mpc_get_prec(x) < wp)
    //     mpc_set_prec(x, wp);

    /* Set 4 * machine precision in u */
    rdpe_set_2dl(u, 1.0, 2 - wp);

    /*    if (MPS_DENSITY_IS_SPARSE(s->active_poly->density)) {  //this line
    fails if enabled? temp_print_sparse(s, p);
        }
    else*/
    std::size_t local_prec =
        std::max(mpf_get_prec(p->mfpc[0]->r), mpf_get_prec(p->mfpc[0]->i));
    local_prec = std::max(local_prec, wp);
    mpfr_set_default_prec(local_prec);
    mpfr_t mpfrre;
    mpfr_t mpfrim;
    arf_t arfre;
    arf_t arfim;
    mpfr_init2(mpfrre, local_prec);
    mpfr_init2(mpfrim, local_prec);
    arf_init(arfre);
    arf_init(arfim);
    ACB acb_coeff(0, 0, local_prec);
    ARB arbre(0, local_prec), arbim(0, local_prec);

    acb_poly_fit_length(pol, MPS_POLYNOMIAL(p)->degree + 1);

    for (int j = 0; j <= MPS_POLYNOMIAL(p)->degree; j++) {
      mpfr_set_f(mpfrre, (p->mfpc[j])->r, MPFR_RNDN);
      mpfr_set_f(mpfrim, (p->mfpc[j])->i, MPFR_RNDN);
      arf_set_mpfr(arfre, mpfrre);
      arf_set_mpfr(arfim, mpfrim);
      arb_set_arf(arbre.r, arfre);
      arb_set_arf(arbim.r, arfim);
      acb_set_arb_arb(acb_coeff.c, arbre.r, arbim.r);
      acb_poly_set_coeff_acb(pol, j, acb_coeff.c);

      // mpc_out_str_2(stdout, 10, 0, 0, p->mfpc[j]);
      // std::cout << std::endl;
    }
    mpfr_clear(mpfrre);
    mpfr_clear(mpfrim);
    arf_clear(arfre);
    arf_clear(arfim);
    /*{
        // mps_with_lock(p->mfpc_mutex[MPS_POLYNOMIAL(p)->degree],
        //             mpc_set(value, p->mfpc[MPS_POLYNOMIAL(p)->degree]););
        mpc_out_str_2(stdout, 10, 0, 0, p->mfpc[MPS_POLYNOMIAL(p)->degree]);
        std::cout << std::endl;

        for (j = MPS_POLYNOMIAL(p)->degree - 1; j >= 0; j--) {
            // mpc_mul_eq(value, x);

            pthread_mutex_lock(&p->mfpc_mutex[j]);
            //            mpc_add_eq(value, p->mfpc[j]);
            std::cout << std::endl;
            pthread_mutex_unlock(&p->mfpc_mutex[j]);
        }
    }*/
  }
  ~ComplexPoly() {
    if (intprec > 0) {
      acb_poly_clear(pol);
    }
  }
  void addRoot(const ACB &c) {
    acb_poly_t temp;
    acb_poly_init(temp);
    acb_poly_set_coeff_si(temp, 1, -1);
    acb_poly_set_coeff_acb(temp, 0, c.c);
    acb_poly_neg(temp, temp);
    if (acb_poly_is_zero(pol) != 0) {
      acb_poly_set(pol, temp);
    } else {
      acb_poly_mul(pol, pol, temp, intprec);
    }
    acb_poly_clear(temp);
  }
  void print(int digits = 5) const {
    acb_poly_printd(pol, digits);
    std::cout << std::endl;
  }
  ComplexPoly graeffe() {
    ComplexPoly temp(intprec);
    acb_poly_graeffe_transform(temp.pol, pol, intprec);
    return temp;
  }
  ComplexPoly derivative() {
    ComplexPoly temp(intprec);
    acb_poly_derivative(temp.pol, pol, intprec);
    return temp;
  }
  void unsquare() { // replaces even exponents x^(2n) with x^n, throws out odd
                    // exponents
    slong len = acb_poly_length(pol);
    if (len == 0) {
      return;
    }
    slong deg = len - 1;
    assert((deg) % 2 == 0);
    ACB tcoeff(0, 0, intprec);
    for (slong i = 0; i <= deg; i += 2) {
      acb_poly_get_coeff_acb(tcoeff.c, pol, i);
      acb_poly_set_coeff_acb(pol, i / 2, tcoeff.c);
    }
    acb_poly_truncate(pol, deg / 2 + 1);
  }
  void explicitCopy(const ComplexPoly &a) {
    acb_poly_set(pol, a.pol);
    intprec = a.intprec;
  }

  ACB getCoeff(slong k) const {
    ACB res(0, 0, intprec);
    acb_poly_get_coeff_acb(res.c, this->pol, k);
    return res;
  }

  void getCoeffMPC(mpc_t &GMPComplex, slong k, slong prec) const {
    ACB coeff = this->getCoeff(k);
    auto rearf = arb_midref(acb_realref(coeff.c));
    auto imarf = arb_midref(acb_imagref(coeff.c));
    mpfr_t rempfr;
    mpfr_t immpfr;
    mpfr_init2(rempfr, prec);
    mpfr_init2(immpfr, prec);
    arf_get_mpfr(rempfr, rearf, MPFR_RNDN);
    arf_get_mpfr(immpfr, imarf, MPFR_RNDN);
    mpf_t regmp;
    mpf_t imgmp;
    mpf_init2(regmp, prec);
    mpf_init2(imgmp, prec);
    mpfr_get_f(regmp, rempfr, MPFR_RNDN);
    mpfr_get_f(imgmp, immpfr, MPFR_RNDN);
    mpc_set_f(GMPComplex, regmp, imgmp);
    mpfr_clear(rempfr);
    mpfr_clear(immpfr);
    mpf_clear(regmp);
    mpf_clear(imgmp);
  }
  slong degree() const { return acb_poly_degree(this->pol); }

  ComplexPoly operator+=(const ComplexPoly &other) {
    acb_poly_add(this->pol, this->pol, other.pol, intprec);
    return *this;
  }
  ComplexPoly operator-=(const ComplexPoly &other) {
    *this = *this - other;
    return *this;
  }
  /*ACB square() {
      acb_sqr(this->c, this->c, intprec);
      return *this;
  }*/
  ComplexPoly operator*=(const ComplexPoly &other) {
    acb_poly_mul(this->pol, this->pol, other.pol, intprec);
    return *this;
  }
  ComplexPoly operator*(const ComplexPoly &other) {
    return ComplexPoly(*this) *= other;
  }
  ComplexPoly operator+(const ComplexPoly &other) {
    return ComplexPoly(*this) += other;
  }
  ComplexPoly operator-(const ComplexPoly &other) {
    ComplexPoly temp(*this);
    ComplexPoly rhs(other); // TODO(orebas) unecessary copy.
    acb_poly_neg(rhs.pol, rhs.pol);
    temp += rhs;
    return temp;
  }

  // ComplexPoly from_mps_poly(mps_polynomial *local_poly) {
  // }
  std::vector<ACB> ArbLibSolve(slong prec) const {
    acb_ptr roots;

    const slong n = this->degree();
    roots = _acb_vec_init(n);
    acb_poly_find_roots(roots, pol, nullptr, 0, prec);
    std::vector<ACB> newroots(n, ACB(0.0, 0.0, prec));
    // could potentially avoid a copy here by using .data() of a vector
    for (slong i = 0; i < n; i++) {
      acb_set(newroots[i].c, roots + i);
    }
    _acb_vec_clear(roots, n);

    sortRootVector(newroots);
    return newroots;
  }

  std::vector<ACB> MPSolve(slong prec) const {
    slong n = this->degree();
    mps_context *status = mps_context_new();
    mps_monomial_poly *poly = mps_monomial_poly_new(status, n);
    mps_context_set_input_prec(status, prec + 1);

    poly->methods.density = MPS_DENSITY_DENSE;
    mps_monomial_poly_raise_precision(status, MPS_POLYNOMIAL(poly), prec);
    //  Set the coefficients. We will solve x^n - 1 in here
    //  mps_monomial_poly_set_coefficient_int (status, poly, 0, -1, 0);
    //  mps_monomial_poly_set_coefficient_int (status, poly, n, 1, 0);
    // std::cout << "monpoly prec" << mps_monomial_poly_get_precision(status,
    // poly)
    //        << std::endl;
    for (slong i = 0; i <= n; i++) {
      mpc_t GMPComplex;
      mpc_init2(GMPComplex, prec);
      this->getCoeffMPC(GMPComplex, i, prec);
      mps_monomial_poly_set_coefficient_f(status, poly, i, GMPComplex);
      mpc_clear(GMPComplex);
    }
    ComplexPoly testc(status, poly, prec);
    // testc.print(50);
    // std::cout << "The diff is: ";
    //((*this) - testc).print(50);
    // std::cout << std::endl;
    //  Select some common output options, i.e. 512 bits of precision
    //  (more or less 200 digits guaranteed) and approximation goal.
    //  Solve the polynomial
    mps_context_set_input_poly(status, MPS_POLYNOMIAL(poly));
    mps_polynomial_set_input_prec(status, MPS_POLYNOMIAL(poly),
                                  prec * 2); // TODO MAGIC NUMBER
    mps_context_set_output_prec(status,
                                prec *
                                    2); // TODO fix decimal vs binary everyhwere
    // MPS_ALGORITHM_SECULAR_GA
    mps_context_select_algorithm(status, MPS_ALGORITHM_STANDARD_MPSOLVE);
    mps_phase phase = no_phase;
    mps_context_set_starting_phase(status, phase);
    mps_context_set_output_goal(status, MPS_OUTPUT_GOAL_APPROXIMATE);

    mps_mpsolve(status);
    // Get the roots in a <code>cplx_t</code> vector. Please note that
    // this make completely useless to have asked 512 bits of output
    // precision, and you should use mps_context_get_roots_m() to get
    // multiprecision approximation of the roots.
    mpc_t *results = nullptr;
    mps_context_get_roots_m(status, &results,
                            nullptr); // inclusion radii discarded.
    // Free the data used. This will free the monomial_poly if you have
    // not done it by yourself.

    for (slong i = 0; i < n; i++) {
    }

    std::vector<ACB> roots;
    for (slong i = 0; i < n; i++) {
      roots.push_back(ACB(results[i], prec));
      // mpc_clear(results[i]);
    }
    mps_monomial_poly_free(status, MPS_POLYNOMIAL(poly));
    mps_context_free(status);
    mpc_vclear(results, n);

    sortRootVector(roots);
    return roots;
  }
};

template <class T = double>
ComplexPoly polyFromRoots(std::vector<T> vec, slong precision) {
  auto length = vec.size();
  auto ptr = _acb_vec_init(length);
  for (std::size_t i = 0; i < length; i++) {
    ACB cc((T)vec[i], (T)0.0, precision);
    acb_set(ptr + i, cc.c); // this is ugly.
  }
  ComplexPoly p(precision);
  acb_poly_product_roots(p.pol, ptr, length, precision);
  _acb_vec_clear(ptr, length);
  return p;
}

ARB inline L2Norm(std::vector<ACB> v1, std::vector<ACB> v2, slong prec) {
  std::size_t n = v1.size();
  assert(v2.size() == n);
  if (n == 0) {
    return ARB(0, prec);
  }
  ARB result = ARB(0.0, prec);
  for (std::size_t i = 0; i < n; i++) {
    ARB temp = (v1[i] - v2[i]).abs();
    result += temp * temp;
  }
  return result;
}
