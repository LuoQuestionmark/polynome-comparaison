#include <Eigen/Eigen>
#include <bits/types/FILE.h>
#include <cassert>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <libcmaes/cmaes.h>
#include <libcmaes/cmastopcriteria.h>
#include <libcmaes/esoptimizer.h>
#include <limits>
#include <random>

#include "evaluation.h"
#include "regulator.h"

using namespace libcmaes;
using namespace Eigen;

// Source - https://stackoverflow.com/a/2704552
// Posted by rep_movsd, modified by community. See post 'Timeline' for change
// history Retrieved 2026-06-07, License - CC BY-SA 4.0
double fRand(double fMin, double fMax) {
    double f = (double)rand() / RAND_MAX;
    return fMin + f * (fMax - fMin);
}

FitFunc fit = [](const double *x, const int N) {
    assert(N == 13);

    MatA_t A;
    VecB_t b;
    Vector<double, 11> y;
    for (int i = 0; i < 11; i++) {
        y[i] = x[i + 2];
    }

    generate_regulator(A, b, x[0], x[1]);
    // generate_regulator(A, b, -1, -1);
    double val = (A * y - b).norm();

    // std::cout << A << std::endl;
    // std::cout << b << std::endl;

    return val;
};

bool launch_solver(int lambda, int max_iteration) {
    std::default_random_engine
        random_engine; // for init variable x with random value

    const int dim = 13;               // size of x
    std::vector<double> x0(dim, 0.0); // declaration of x

    // random value init
    x0[0] = fRand(-1, 0);
    x0[1] = fRand(-10, 0);
    for (int i = 2; i < 10; i++) {
        x0[i] = fRand(0, 1);
    }
    for (int i = 10; i < 13; i++) {
        x0[i] = fRand(-5, 5);
    }

    double sigma = 0.1;

    // init bounds
    double lbounds[dim], ubounds[dim];
    for (int i = 0; i < 2; i++) {
        lbounds[i] = -std::numeric_limits<double>::infinity();
        ubounds[i] = 0;
    }
    for (int i = 2; i < 10; i++) {
        lbounds[i] = 0;
        ubounds[i] = std::numeric_limits<double>::infinity();
    }
    for (int i = 10; i < 13; i++) {
        lbounds[i] = -std::numeric_limits<double>::infinity();
        ubounds[i] = std::numeric_limits<double>::infinity();
    }

    GenoPheno<pwqBoundStrategy> gp(lbounds, ubounds, dim); // bounds

    // int lambda = 100; // offsprings at each generation.
    // CMAParameters<> cmaparams(x0, sigma);
    CMAParameters<GenoPheno<pwqBoundStrategy>> cmaparams(dim, &x0.front(),
                                                         sigma, lambda, 0, gp);
    // cmaparams.set_max_iter(max_iteration);
    // cmaparams.set_stopping_criteria(CMAStopCritType::TOLHISTFUN, false);
    // cmaparams.set_stopping_criteria(CMAStopCritType::STAGNATION, false);
    // cmaparams.set_stopping_criteria(CMAStopCritType::EQUALFUNVALS, false);
    // cmaparams.set_stopping_criteria(CMAStopCritType::, false);
    cmaparams.set_ftarget(1e-6);
    cmaparams.set_quiet(false);
    cmaparams.set_algo(IPOP_CMAES);

    // cmaparams._algo = BIPOP_CMAES;
    CMASolutions cmasols = cmaes<GenoPheno<pwqBoundStrategy>>(fit, cmaparams);

    std::cout << "best solution\n"
              << gp.pheno(cmasols.get_best_seen_candidate().get_x_dvec())
              << std::endl;
    std::cout << "optimization took " << cmasols.elapsed_time() / 1000.0
              << " seconds\n";
    std::cout << cmasols.status_msg() << std::endl;

    cmasols.sort_candidates();

    auto candidate = cmasols.candidates()[0].get_x_pheno_dvec(cmaparams);
    fit(candidate.data(), 13);

    double best = cmasols.get_candidate(0).get_fvalue();
    double worst =
        cmasols.get_candidate(cmasols.candidates().size() - 1).get_fvalue();

    printf("%.8f, %.8f\n", best, worst);

    Vector<double, 13> best_vector(
        gp.pheno(cmasols.get_best_seen_candidate().get_x_dvec()));

    if (tau_eval(best_vector, 1)) {
        puts("stable");
    } else {
        puts("not stable");
    }

    return cmasols.run_status() == 1;
}

int main(int argc, char **argv) {
    const int total_test = 1;
    int total_success    = 0;

    int population_aka_lambda = 500;
    int max_iteration         = 1000;

    for (int i = 0; i < total_test; i++) {
        if (launch_solver(population_aka_lambda, max_iteration) == true) {
            total_success += 1;
        }
    }

    printf("success rate: %f\n", (float)total_success / (float)total_test);
}
