

#include <chrono>
#include <iostream>
#include <cmath>

#include <functional>
#include <memory>
#include <utility>
#include <filesystem>

#include "arb.h"
#include "acb_poly.h"
#include "acb_calc.h"

#include "acb.h"
#include "arb_fmpz_poly.h"
#include "flint/arith.h"
#include "flint/profiler.h"

#define _MPS_PRIVATE
#include <mps/mps.h>
#include <string.h>

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#define MPSOLVE_GETOPT_STRING "a:G:D:d::t:o:O:j:S:O:i:vl:bp:rs:c"

#ifndef __WINDOWS
#include <signal.h>

void status(int signal, mps_context *s);

#undef _MPS_PRIVATE
#endif

#include <benchmark.hpp>
//#include "Hungarian.h"
// void * cleanup_context(mps_context *ctx, void *user_data);

int main(int argc, char **argv)
{
    // swap debugging
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
        // std::cout << " " << a.intprec << " " << b.intprec << " " << temp.intprec << std::endl;

        // std::cout << a << " " << b << " " << temp;
        // std::cout << " " << a.intprec << " " << b.intprec << " " << temp.intprec << std::endl;

        //    std::cout << " " << a << " " << b;
        //   std::cout << " " << a.intprec << " " << b.intprec << " " << std::endl;
    */
    // end swap debugging

    std::vector<int> dvec;
    for (int i = 1; i <= 40; i++)
    {
        dvec.push_back(i);
    }
    auto p = polyFromRoots<int>(dvec, 800);
    p.print(25);
    auto res = p.MPSolve(800);
    std::cout << "Roots:\n"
              << res << std::endl
              << std::endl;
    return 0;
    for (auto const &dir_entry : std::filesystem::directory_iterator{"."})
    {
        if (dir_entry.path().extension() == ".pol")
        {
            // std::cout << dir_entry << " " << dir_entry.path().extension() << std::endl;
            parsePol(dir_entry.path().string());
        }
        // std::cout << dir_entry << " " << dir_entry.path().extension() << std::endl;
    }
}

void *cleanup_context(mpscp ctx, void *user_data, mps_polynomial *poly, mpscp &s, bool outp)
{
    /* Check for errors */
    if (mps_context_has_errors(ctx))
    {
        mps_print_errors(ctx);
        return NULL;
    }

    /* Output the roots */
    if (outp)
    {
        mps_output(ctx);
    }
    /* Free used data */
    if (poly)
        mps_polynomial_free(ctx, poly);
    mps_context_free(ctx);

    s = NULL;

    return NULL;
}

void parsePol(const std::string &polfilename)
{
    mps_context *local_s = NULL;
    mps_polynomial *local_poly = NULL;

    /* Create a new status */
    local_s = mps_context_new();

    /* Associate the SIGUSR1 signal to the mps_dump () function */
    //#ifndef __WINDOWS
    //    signal(SIGUSR1, status);
    //#endif

    /* Leave this value to -1 if not precision has been enforced
     * by the user. Otherwise, override the polynomial input
     * precision. */
    long int input_precision = -1;
    char *inline_poly = NULL;

    FILE *infile;

    /* Parse options */
    // mps_opt *opt = NULL; //commented out but we might need it later
    mps_phase phase = no_phase;

    /* This will be true if the user has explicitely selected the algorithm,
       otherwise we will be using our own heuristic. */
    mps_boolean explicit_algorithm_selection = false;
    // opt = NULL;
    /* we could set a logfile with the following code:
   {
                FILE *logstr = fopen(opt->optvalue, "w");
                if (!logstr)
                    mps_error(s, "Cannot open selected log file.");
                mps_context_set_log_stream(s, logstr);
            }
     */
    // we could set output format with mps_context_set_output_format //see main() for examples
    infile = fopen(polfilename.c_str(), "r");
    if (!infile)
    {
        mps_error(local_s, "Cannot open input file for read, aborting.");
        mps_print_errors(local_s);
        return; // EXIT_FAILURE;
    }

    /* Parse the input stream and if a polynomial is given as output,
     * allocate also a secular equation to be used in regeneration */
    local_poly = mps_parse_stream(local_s, infile);

    if (!local_poly)
    {
        mps_error(local_s, "Error while parsing the polynomial, aborting.");
        mps_print_errors(local_s);
        return; // EXIT_FAILURE;
    }
    else
        mps_context_set_input_poly(local_s, local_poly);
    /*local_s->active_poly = local_poly;
    std::cout << local_s << std::endl;
    std::cout << local_s->active_poly << std::endl;
    std::cout << local_poly << std::endl;
*/
    /* Override input precision if needed */
    if (input_precision >= 0)
        mps_polynomial_set_input_prec(local_s, local_poly, 2000); // TODO MAGIC NUMBER

    /* Perform some heuristic for the algorithm selection, but only if the user
     * hasn't explicitely selected one. */
    if (!explicit_algorithm_selection)
    {
        mps_context_select_algorithm(local_s, (MPS_IS_MONOMIAL_POLY(local_poly) &&
                                               MPS_DENSITY_IS_SPARSE(local_poly->density))
                                                  ? MPS_ALGORITHM_STANDARD_MPSOLVE
                                                  : MPS_ALGORITHM_SECULAR_GA);
    }

    /* Close the file if it's not stdin */
    if (!inline_poly)
        fclose(infile);

    /* Select the starting phase according to user input */
    mps_context_set_starting_phase(local_s, phase);

    // BEGIN WILK20 DEBUGGING CODE
    /*
    mps_monomial_poly *local_poly2 = MPS_POLYNOMIAL_CAST(mps_monomial_poly, local_poly);

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

        arb_fmpz_poly_complex_roots(roots, fac->p + i, flags, compd * 3.32193 + 2);

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

void temp_print(mps_context *local_s, mps_polynomial *poly)
{
    int i;
    // mps_monomial_poly *p = MPS_MONOMIAL_POLY(poly);
    mps_monomial_poly *p = (mps_monomial_poly *)poly;

    std::cout << poly->degree << std::endl;
    fflush(stdout);

    for (i = 0; i <= poly->degree; i++)
    {
        rdpe_out(p->dap[i]);
        std::cout << p->fap[i] << std::endl;
        std::cout << p->fpr[i] << std::endl;
        cplx_out(p->fpc[i]);
        rdpe_out(p->dpr[i]);
        cdpe_out(p->dpc[i]);
        mpf_out_str(stdout, 10, 0, p->mfpr[i]);
        mpc_out_str(stdout, 10, 0, p->mfpc[i]);

        if (i < poly->degree)
        {
            mpc_out_str(stdout, 10, 0, p->mfppc[i]);
        }
        mpq_out_str(stdout, 10, p->initial_mqp_r[i]);
        mpq_out_str(stdout, 10, p->initial_mqp_i[i]);
        std::cout << p->spar[i] << std::endl;
    }
}

std::vector<ACB> powersums(slong depth, slong k, slong prec)
{
    std::vector<ACB> results;
    int i, d;
    for (d = 0; d < depth + 1; d++)
    {
        ACB temp(0, 0, prec);

        for (i = 1; i <= k; i++)
        {
            temp += (ACB((int)i, (int)0, prec)).power(d);
        }
        results.push_back(temp);
    }
    return results;
}
void solveCompare(mps_context *local_s, mps_polynomial *poly_local_poly, const std::string &polfilename)
{
    mps_monomial_poly *local_poly = MPS_POLYNOMIAL_CAST(mps_monomial_poly, poly_local_poly);

    std::cout << "Starting " << polfilename << std::endl;
    mps_context_set_output_prec(local_s, 512);
    mps_context_set_output_goal(local_s, MPS_OUTPUT_GOAL_APPROXIMATE);
    /* Solve the polynomial */
    auto funcMPS = [&]() -> void
    {
        mps_mpsolve(local_s);
    };

    // ComplexPoly new_poly;
    //  auto parsed_poly = ComplexPoly::from_mps_poly(mps_polynomial * local_poly);

    slong prec = 2000; // TODO MAGIC NUMBER
    ComplexPoly acb_style_poly = ComplexPoly(local_s, local_poly, prec);
    // std::cout << "The poly might be ";
    // acb_style_poly.print();
    std::cout << std::endl;
    std::vector<ComplexPoly> pvec, qvec;
    // DLGFG(acb_style_poly, pvec, qvec, 90);
    int depth = 20;
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
    auto funcRS1 = [&]() -> void
    {
        RootSolver1Roots = RS1.solve();
    };
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
    // uses recursive FG and DLG to estimate largest and smallest roots, and root radii

    double mpsTime = measure<std::chrono::milliseconds>::execution(funcMPS); // runs mpsolve
    double RS1Time = measure<std::chrono::milliseconds>::execution(funcRS1); // runs mpsolve

    fflush(stdout);
    mpc_t *results = NULL;
    // std::cout << "MPS SOLVE" << std::endl;
    mps_context_get_roots_m(local_s, &results, NULL);
    std::vector<ACB> mps_roots;
    for (slong rt = 0; rt < acb_style_poly.degree(); rt++)
    {
        mps_roots.push_back(ACB(results[rt], prec));
    }

    // std::cout << "MPS SOLVE1" << std::endl;
    //  std::cout << mps_roots;
    struct
    {
        bool operator()(const ACB &a, const ACB &b) const
        {
            // ARF l = a.abs_ubound_arf();
            // ARF r = b.abs_ubound_arf();
            // int i = arf_cmp(l.f, r.f);
            // return i < 0;
            if (arb_lt(a.abs().r, b.abs().r))
            {
                return true;
            }
            if (arb_gt(a.abs().r, b.abs().r))
            {
                return false;
            }
            ARB arg1(0.0, a.intprec);
            ARB arg2(0.0, b.intprec);
            acb_arg(arg1.r, a.c, a.intprec);
            acb_arg(arg2.r, b.c, a.intprec);
            if (arb_lt(arg1.r, arg2.r))
            {
                return true;
            }
            if (arb_gt(arg1.r, arg2.r))
            {
                return false;
            }
            return false;
        }
    } customLess;

    struct
    {
        bool operator()(const ACB &a, const ACB &b) const
        {
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

            if (arb_lt(re1.r, re2.r))
            {
                return true;
            }
            if (arb_gt(re1.r, re2.r))
            {
                return false;
            }
            if (arb_lt(im1.r, im2.r))
            {
                return true;
            }
            if (arb_gt(im1.r, im2.r))
            {
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
            costMatrix[i].push_back((mps_roots[i] - RootSolver1Roots[j]).abs());
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

    std::cout
        << "RootSolver1" << std::endl;
    std::cout << RootSolver1Roots;

    std::cout << "MPS SOLVE2" << std::endl;
    std::cout << mps_roots;

    std::cout << polfilename << " took " << mpsTime << "ms. for MPS " << std::endl;
    std::cout << polfilename << " took " << RS1Time << "ms. for RS1." << std::endl;
    // auto res = ACBVectorComp(mps_roots, rootapprox.back());
    // std::cout << res << std::endl;
    std::cout << mps_roots.size() << std::endl;

    auto rootapprox = RS1.rootapprox;

    std::cout << rootapprox.back().size() << std::endl;

    for (int ri = 1; ri < rootapprox.size(); ri++)
    {
        std::sort(rootapprox[ri].begin(), rootapprox[ri].end(), customLessLex);
        std::cout
            << "\nRow " << ri << " " << ACBVectorComp(mps_roots, rootapprox[ri]);
    }
    mpc_vfree(results);
    cleanup_context(local_s, NULL, MPS_POLYNOMIAL(local_poly), local_s, false);
}

void solveDebug(mps_context *local_s, mps_polynomial *poly_local_poly, const std::string &polfilename)
{
    mps_monomial_poly *local_poly = MPS_POLYNOMIAL_CAST(mps_monomial_poly, poly_local_poly);
    /* Solve the polynomial */
    auto func = [&]() -> void
    {
        mps_mpsolve(local_s);
    };

    // ComplexPoly new_poly;
    //  auto parsed_poly = ComplexPoly::from_mps_poly(mps_polynomial * local_poly);
    slong prec = 2000; // TODO MAGIC NUMBER
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

    for (int depth = 0; depth < 15; depth++)
    {
        std::cout << depth << "[min]" << RRE.minRadius(depth);
        std::cout << "[max]" << RRE.maxRadius(depth);
    }

    for (int depth = 0; depth < 15; depth++)
    {
        std::cout << depth << "[minRoot]" << RRE.minRoot(depth);

        std::cout << depth << "[minRoot]" << RRE.maxRoot(depth);
    }

    double initTime = measure<std::chrono::milliseconds>::execution(func);
    fflush(stdout);
    std::cout << polfilename << " took " << initTime << "ms." << std::endl;
    cleanup_context(local_s, NULL, MPS_POLYNOMIAL(local_poly), local_s, false);
}

class integrand_parameters
{
public:
    ComplexPoly p;
    ARB r;
    slong xpower;
    integrand_parameters(const ComplexPoly &initpol, ARB radius, slong pow) : p(initpol.intprec), r(0, radius.intprec)
    {
        p.explicitCopy(initpol);
        arb_set(r.r, radius.r);
        xpower = pow;
    }
};

int pprime_over_p_integrand(acb_ptr out, const acb_t inp, void *paramsv, slong order, slong prec)
{
    // this is the integrand to evaluate integral p'/p, order-th derivative.  inp goes from 0 to 1.

    ACB a(0, 0, prec), b(0, 0, prec), z(0, 0, prec);
    ACB inp_copy(0, 0, prec);
    integrand_parameters *params = reinterpret_cast<integrand_parameters *>(paramsv);
    acb_set(inp_copy.c, inp);
    acb_mul_2exp_si(inp_copy.c, inp_copy.c, 1);      // inp *=2;
    acb_exp_pi_i(z.c, inp_copy.c, inp_copy.intprec); // z = exp(2 * pi *i * inp)
    acb_mul_arb(z.c, z.c, params->r.r, z.intprec);   // z = r e ^(2 * pi * inp);

    acb_poly_evaluate2(a.c, b.c, params->p.pol, z.c, params->p.intprec);
    acb_div(out, b.c, a.c, b.intprec); // out = p'/p at re^(2 pi i inp)
    acb_pow_si(z.c, z.c, params->xpower, z.intprec);
    acb_div(out, out, z.c, z.intprec); // return p'(z)/p(z)z^n

    return 0;
}

// below approximates first depth derivatives of pprime/p via cauchy integral formula
void DerivEst(std::vector<ACB> &results, const ComplexPoly &p, ARB minrad, int depth)
{
    mag_t tol;
    slong prec, goal;
    acb_calc_integrate_opt_t options;

    acb_calc_integrate_opt_init(options);

    prec = p.intprec; // p.intprec;
    goal = prec;
    mag_init(tol);
    mag_set_ui_2exp_si(tol, 1, -prec);

    long k = 0;
    ACB deriv_est(0, 0, prec);
    results.resize(depth + 1, ACB(0, 0, prec));
    integrand_parameters params(p, minrad, k);
    ACB start(0, 0, prec), end(1, 0, prec);
    auto func(&pprime_over_p_integrand);
    for (; k <= depth; k++)
    {
        params.xpower = k;
        acb_calc_integrate(deriv_est.c, func, &params, start.c, end.c, goal, tol, options, prec);
        fmpz_t t;
        fmpz_fac_ui(t, k);
        acb_mul_fmpz(deriv_est.c, deriv_est.c, t, deriv_est.intprec);
        results[k].explicitCopy(deriv_est);
    }
}
