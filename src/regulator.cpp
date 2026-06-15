#include "regulator.h"
#include <Eigen/src/Core/Matrix.h>

void generate_regulator(MatA_t &A, VecB_t &b, double a1, double a2) {
    Matrix2d Ae;
    Vector<double, 3> be;
    Matrix<double, 2, 3> A1, A2, A3, A4;
    Matrix<double, 3, 3> A5;

    Ae << 1, 1, -a1, -a2;
    A1 << 1, -1, 0, 1, 0, -1;
    A2 << 1, +1, 0, 2, 0, +1;
    A3 << 1, -1, 0, 2, 0, -1;
    A4 << 1, +1, 0, 1, 0, +1;
    A5 << 0, 2, 0, 0, 1, 2, 0, 0, 1;

    be << 1, -(a1 + a2), a1 * a2;

    // assignment of A
    A.setZero();

    for (int i = 0; i < 4; i++) {
        A.block(i * 2, i * 2, 2, 2) = Ae;
        A.block(8 + i * 3, 8, 3, 3) = A5;
    }

    A.block(0, 8, 2, 3) = -A1;
    A.block(2, 8, 2, 3) = -A2;
    A.block(4, 8, 2, 3) = -A3;
    A.block(6, 8, 2, 3) = -A4;

    // assignment of b
    b.setZero();

    for (int i = 0; i < 4; i++) {
        b.segment(8 + i * 3, 3) = be;
    }
}
