#ifndef REGULATOR_H
#define REGULATOR_H

#include <Eigen/Eigen>

using namespace Eigen;

typedef Matrix<double, 20, 11> MatA_t;
typedef Vector<double, 20> VecB_t;

void generate_regulator(MatA_t &A, VecB_t &b, double a1, double a2);

#endif
