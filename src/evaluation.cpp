#include "evaluation.h"
#include "regulator.h"
#include <Eigen/src/Core/Matrix.h>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <unsupported/Eigen/Polynomials>

// correctness of this function has been verified
Vector<double, 5> conv_vec_3_3(const Vector<double, 3> &a,
                               const Vector<double, 3> &b) {
    Vector<double, 5> out;
    out.setZero();

    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            out[i + j] += a[i] * b[j];
        }
    }

    return out;
}

bool tau_eval(const Vector<double, 13> &x, double tau) {
    double c_min = 1;
    double den_0 = 1;
    double den_2 = 2;

    double a_min = 1.25 - 0.25 * tau;
    double a_max = 1.25 + 0.25 * tau;

    double b_min = -tau;
    double b_max = +tau;

    Vector<double, 3> N1, N2, N3, N4, D1, D2, D3, D4;
    Vector<double, 3> Nc, Dc;
    Vector<double, 5> P1, P2, P3, P4;

    N1 << c_min, 0, a_min;
    D1 << den_2, b_min, den_0;

    N2 << c_min, 0, a_max;
    D2 << den_2, b_min, den_0;

    N3 << c_min, 0, a_min;
    D3 << den_2, b_max, den_0;

    N4 << c_min, 0, a_max;
    D4 << den_2, b_max, den_0;

    Nc << 0, x[10], 0;
    Dc << x[11], 0, x[12];

    P1 = conv_vec_3_3(N1, Nc) + conv_vec_3_3(D1, Dc);
    P2 = conv_vec_3_3(N2, Nc) + conv_vec_3_3(D2, Dc);
    P3 = conv_vec_3_3(N3, Nc) + conv_vec_3_3(D3, Dc);
    P4 = conv_vec_3_3(N4, Nc) + conv_vec_3_3(D4, Dc);

    // see
    // https://stackoverflow.com/questions/61074378/calculating-roots-of-polynom-in-c-with-eigen-library#61075738
    // polynomial solver functions just fine; only difference is the order of
    // coefficients, it is the reversed order compare to MATLAB impl
    Eigen::PolynomialSolver<double, Eigen::Dynamic> solver;

    solver.compute(P1.reverse());
    if (solver.roots().real().cwiseSign().cwiseMax(0).any()) {
        std::cout << "P1: " << P1 << std::endl;
        return false;
    }

    solver.compute(P2.reverse());
    if (solver.roots().real().cwiseSign().cwiseMax(0).any()) {
        std::cout << "P2: " << P2 << std::endl;
        return false;
    }

    solver.compute(P3.reverse());
    if (solver.roots().real().cwiseSign().cwiseMax(0).any()) {
        std::cout << "P3: " << P3 << std::endl;
        return false;
    }

    solver.compute(P4.reverse());
    if (solver.roots().real().cwiseSign().cwiseMax(0).any()) {
        std::cout << "P4: " << P4 << std::endl;
        return false;
    }

    return true;
}
