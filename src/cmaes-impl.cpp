#include <Eigen/Eigen>
#include <asm-generic/errno-base.h>
#include <bits/types/FILE.h>
#include <cassert>
#include <cerrno>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <libcmaes/cmaes.h>
#include <libcmaes/cmastopcriteria.h>
#include <limits>
#include <random>
#include <vector>

#include "nlohmann/json.hpp" // config using json file

#include "evaluation.h"
#include "regulator.h"

using namespace libcmaes;
using namespace Eigen;
using json = nlohmann::json;

FILE *GLOBAL_LOG_FILE_PTR = NULL;

// Source - https://stackoverflow.com/a/2704552
// Posted by rep_movsd, modified by community. See post 'Timeline' for change
// history Retrieved 2026-06-07, License - CC BY-SA 4.0
double fRand(double fMin, double fMax) {
    double f = (double)rand() / RAND_MAX;
    return fMin + f * (fMax - fMin);
}

// fitness function, using the norm2 value of $Ay - b$ as the goal
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

// action of launching the solver a single time;
// this function will be called multiple times for one configuration to measure
// the average performance
bool launch_solver(int lambda, int max_iteration, double fitness_target) {
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
    cmaparams.set_ftarget(fitness_target);
    // cmaparams.set_quiet(false);
    // cmaparams.set_algo(IPOP_CMAES);
    // cmaparams._algo = BIPOP_CMAES;

    CMASolutions cmasols = cmaes<GenoPheno<pwqBoundStrategy>>(fit, cmaparams);

    // std::cout << "best solution\n"
    //           << gp.pheno(cmasols.get_best_seen_candidate().get_x_dvec())
    //           << std::endl;
    // std::cout << "optimization took " << cmasols.elapsed_time() / 1000.0
    //           << " seconds\n";
    // std::cout << cmasols.status_msg() << std::endl;

    cmasols.sort_candidates();

    double best = cmasols.get_candidate(0).get_fvalue();
    double worst =
        cmasols.get_candidate(cmasols.candidates().size() - 1).get_fvalue();

    if (GLOBAL_LOG_FILE_PTR) {
        fprintf(GLOBAL_LOG_FILE_PTR, "min err: %.8f, max err: %.8f\n", best,
                worst);
    }

    Vector<double, 13> best_vector(
        gp.pheno(cmasols.get_best_seen_candidate().get_x_dvec()));

    bool stable = tau_eval(best_vector, 1);
    return stable;
}

int main(int argc, char **argv) {
    if (argc <= 1) {
        errno = EINVAL;
        perror("usage: cmaes-impl <config file.json>");
        exit(EXIT_FAILURE);
    }

    printf("Load config file %s...\n", argv[1]);

    std::ifstream config(argv[1]);
    if (!config.good()) {
        errno = ENOENT;
        perror("cannot open json config file");
        exit(EXIT_FAILURE);
    }

    json config_data = json::parse(config);

    const std::vector<std::string> fields = { "total test", "population",
                                              "max iteration", "fitness target",
                                              "log file" };
    for (const std::string &key : fields) {
        if (!config_data.contains(key)) {
            errno = ENOKEY;
            fprintf(stderr, "config file does not contain key: %s\n",
                    key.c_str());
            exit(EXIT_FAILURE);
        }
    }

    int total_test    = config_data["total test"];
    int total_success = 0;

    int population        = config_data["population"];
    int max_iteration     = config_data["max iteration"];
    double fitness_target = config_data["fitness target"];

    puts("Config load OK");

    std::string log     = config_data["log file"];
    GLOBAL_LOG_FILE_PTR = fopen(log.c_str(), "w+");

    if (!GLOBAL_LOG_FILE_PTR) {
        perror("cannot create log file");
        exit(EXIT_FAILURE);
    }

    puts("Log file creation OK");

    if (GLOBAL_LOG_FILE_PTR) {
        fprintf(GLOBAL_LOG_FILE_PTR,
                "total test count: %d, population size: %d, maximum iteration: "
                "%d, fitness target: %.3e\n",
                total_test, population, max_iteration, fitness_target);
    }

    puts("Start testing");

    for (int i = 0; i < total_test; i++) {
        printf("iteration no. %d\n", i + 1);
        if (launch_solver(population, max_iteration, fitness_target) == true) {
            total_success += 1;
        }
    }

    puts("Finish testing");

    if (GLOBAL_LOG_FILE_PTR) {
        fprintf(GLOBAL_LOG_FILE_PTR, "success rate: %f\n",
                (float)total_success / (float)total_test);
        fclose(GLOBAL_LOG_FILE_PTR);
    }
}
