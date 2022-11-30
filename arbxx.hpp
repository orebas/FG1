

#include "arb.h"
#include "acb_poly.h"
#include "acb_calc.h"
#include <cassert>
template <class T = double>
class Array3d {
   private:
    int w;
    int h;
    int d;

   public:
    std::vector<T> data;

    Array3d(int width, int height, int depth) : w(width), h(height), d(depth), data(w * h * d, 0) {}

    inline T &at(int x, int y, int z) { return data[x * h * d + y * d + z]; }

    inline T at(int x, int y, int z) const {
        return data[x * h * d + y * d + z];
    }

    inline T &operator()(int x, int y, int z) {
        return data[x * h * d + y * d + z];
    }

    inline T operator()(int x, int y, int z) const {
        return data[x * h * d + y * d + z];
    }
    inline int width() { return w; }
    inline int height() { return h; }
    inline int depth() { return d; }
};

template <class T = double>
class Array4d {
   private:
    int w;
    int h;
    int d;
    int l;

   public:
    std::vector<T> data;

    Array4d(int r, int s, int t, int u) : w(r), h(s), d(t), l(u), data(w * h * d * l, 0) {}

    inline T &at(int x, int y, int z, int t) {
        return data[x * h * d * l + y * d * l + z * l + t];
    }

    inline T at(int x, int y, int z, int t) const {
        return data[x * h * d * l + y * d * l + z * l + t];
    }

    inline T &operator()(int x, int y, int z, int t) {
        return data[x * h * d * l + y * d * l + z * l + t];
    }

    inline T operator()(int x, int y, int z, int t) const {
        return data[x * h * d * l + y * d * l + z * l + t];
    }

    inline int internalRef(int x, int y, int z, int t) {
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
    ARF(ARF &&other)
    noexcept : intprec(other.intprec) {
        arf_swap(arf_ptr_int(), other.arf_ptr_int());
        other.intprec = -1;
    }

    ARF &operator=(const ARF &other) {
        return *this = ARF(other);
    }

    ARF &operator=(ARF &&other) noexcept {
        // std::swap(c, other.c);
        arf_swap(arf_ptr_int(), other.arf_ptr_int());  // Do I need to deinitialize here?
        return *this;
    }

    ~ARF() {
        if (intprec > 0)
            arf_clear(f);
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
    ARB(const ARB &other) : ARB(other.r, other.intprec) {}
    ARB(ARB &&other)
    noexcept : intprec(other.intprec) {
        arb_swap(arb_ptr_int(), other.arb_ptr_int());
        other.intprec = -1;
    }

    ARB &operator=(const ARB &other) {
        return *this = ARB(other);
    }

    ARB root_ui(slong e) {
        ARB temp(0, intprec);
        arb_root_ui(temp.r, this->r, e, intprec);
        return temp;
    }

    ARB &operator=(ARB &&other) noexcept {
        // std::swap(c, other.c);
        intprec = other.intprec;
        arb_swap(arb_ptr_int(), other.arb_ptr_int());  // Do I need to deinitialize here?
        other.intprec = -1;
        return *this;
    }

    ARB operator/=(const ARB &other) {
        arb_div(this->r, this->r, other.r, intprec);
        return *this;
    }

    ARB operator/(const ARB &other) {
        return ARB(*this) /= other;
    }
    ARB operator-=(const ARB &other) {
        arb_sub(this->r, this->r, other.r, intprec);
        return *this;
    }

    ARB operator-(const ARB &other) {
        return ARB(*this) -= other;
    }

    ARB operator+=(const ARB &other) {
        arb_add(this->r, this->r, other.r, intprec);
        return *this;
    }

    ARB operator+(const ARB &other) {
        return ARB(*this) += other;
    }

    ARB operator*=(const ARB &other) {
        arb_mul(this->r, this->r, other.r, intprec);
        return *this;
    }
    ARB operator*(const ARB &other) {
        return ARB(*this) *= other;
    }

    bool operator<(const ARB &other) {
        return arb_lt(this->r, other.r);
    }
    bool operator>(const ARB &other) {
        return arb_gt(this->r, other.r);
    }

    ~ARB() {
        if (intprec > 0)
            arb_clear(r);
    }
    void explicitCopy(const ARB &a) {
        arb_set(r, a.r);
        intprec = a.intprec;
    }
    void print() const { arb_printd(r, 30); }
};

ARB fabs(ARB &a) {
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

    ACB(int a, int b, slong ip) : intprec(ip) {
        acb_init(c);
        acb_set_si_si(c, a, b);
    }

    ACB(slong a, slong b, slong ip) : intprec(ip) {
        acb_init(c);
        acb_set_si_si(c, a, b);
    }
    ACB(const mpc_t &p, slong local_prec) : intprec(local_prec) {
        acb_init(c);
        mpfr_t mpfrre,
            mpfrim;
        arf_t arfre, arfim;  // TODO Is this a memory leak? add init and clear, maybe.
        mpfr_init2(mpfrre, local_prec);
        mpfr_init2(mpfrim, local_prec);
        arf_init(arfre);
        arf_init(arfim);
        ARB arbre(0, local_prec), arbim(0, local_prec);
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
    ACB(ACB &&other)
    noexcept : intprec(other.intprec) {
        // acb_ptr_int() = 0;
        acb_swap(acb_ptr_int(), other.acb_ptr_int());
        other.intprec = -1;
        // c = other.c;
        // other.c = nullptr;
    }

    ACB &operator=(const ACB &other) {
        return *this = ACB(other);
    }

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
    ACB operator*(const ACB &other) {
        return ACB(*this) *= other;
    }
    ACB operator+(const ACB &other) {
        return ACB(*this) += other;
    }
    ACB operator/(const ACB &other) {
        return ACB(*this) /= other;
    }
    ACB operator-(const ACB &other) {
        return ACB(*this) -= other;
    }
    ACB power(slong p) {
        ACB temp(0, 0, intprec);
        acb_pow_si(temp.c, this->c, p, intprec);
        return temp;
    }
    ARB abs() const {
        assert(this->intprec > 0);
        ARB temp(0, intprec);
        // std::cout << "temp init:" << arb_get_str(temp.r, 20, 0) << std::endl;
        // std::cout << "cplx init:" << arb_get_str(acb_realref(this->c), 20, 0) << " " << arb_get_str(acb_imagref(this->c), 20, 0);
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
        if (intprec > 0)
            acb_clear(c);
    }

    void explicitCopy(const ACB &a) { acb_set(c, a.c); }
    void print(int digits = 5) { acb_printd(c, digits); }
};

// ACB operator*(const ACB &lhs, const ACB &rhs) {
//     return ACB(lhs) *= rhs;
// }

class ComplexPoly {
   public:
    acb_poly_t pol;
    slong intprec = 1;

    typedef acb_poly_struct *acb_poly_ptr;

    acb_poly_ptr acb_poly_ptr_int() { return pol; }

    ComplexPoly(slong ip) : intprec(ip) { acb_poly_init(pol); }
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
    ComplexPoly(mps_context *s, mps_monomial_poly *p, slong wp) : intprec(wp) {
        acb_poly_init(pol);
        // int i;
        // rdpe_t apol, ax
        rdpe_t u;
        // cdpe_t cx;

        pthread_mutex_lock(&p->mfpc_mutex[0]);
        if (mpc_get_prec(p->mfpc[0]) < wp) {
            pthread_mutex_unlock(&p->mfpc_mutex[0]);
            mps_monomial_poly_raise_precision(s, MPS_POLYNOMIAL(p), wp);
        } else
            pthread_mutex_unlock(&p->mfpc_mutex[0]);

        // if (mpc_get_prec(x) < wp)
        //     mpc_set_prec(x, wp);

        /* Set 4 * machine precision in u */
        rdpe_set_2dl(u, 1.0, 2 - wp);

        int j;

        /*    if (MPS_DENSITY_IS_SPARSE(s->active_poly->density)) {  //this line fails if enabled?
                temp_print_sparse(s, p);
            }
        else*/
        slong local_prec = std::max(mpf_get_prec(p->mfpc[0]->r), mpf_get_prec(p->mfpc[0]->i));
        local_prec = std::max(local_prec, wp);
        mpfr_set_default_prec(local_prec);
        mpfr_t mpfrre,
            mpfrim;
        arf_t arfre, arfim;
        mpfr_init2(mpfrre, local_prec);
        mpfr_init2(mpfrim, local_prec);
        arf_init(arfre);
        arf_init(arfim);
        ACB acb_coeff(0, 0, local_prec);
        ARB arbre(0, local_prec), arbim(0, local_prec);

        acb_poly_fit_length(pol, MPS_POLYNOMIAL(p)->degree + 1);

        for (j = 0; j <= MPS_POLYNOMIAL(p)->degree; j++) {
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
        if (acb_poly_is_zero(pol)) {
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
    void unsquare() {  // replaces even exponents x^(2n) with x^n, throws out odd exponents
        slong len, deg;
        len = acb_poly_length(pol);
        if (len == 0) return;
        deg = len - 1;
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

    ACB getCoeff(slong k) {
        ACB res(0, 0, intprec);
        acb_poly_get_coeff_acb(res.c, this->pol, k);
        return res;
    }
    slong degree() const {
        return acb_poly_degree(this->pol);
    }

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
        ComplexPoly rhs(other);  // TODO unecessary copy.
        acb_poly_neg(rhs.pol, rhs.pol);
        temp += rhs;
        return temp;
    }

    // ComplexPoly from_mps_poly(mps_polynomial *local_poly) {
    // }
};

template <class T = double>
ComplexPoly polyFromRoots(std::vector<T> vec,slong precision){
    auto length = vec.size();
    auto ptr = _acb_vec_init(length);
    slong i=0;
    for(;i< length; i++){
        ACB cc(vec[i],0.0,precision);
        acb_set(ptr+i, cc.c); //this is ugly.
    }
    ComplexPoly p(precision);
    acb_poly_product_roots(p.pol,ptr,length,precision);
    _acb_vec_clear(ptr,length);
return p;
}
