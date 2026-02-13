/// @file unit.cpp
/// @brief Source file for unit function
/// @date 2025-04-15
/// @author Martín Hernández Tonzán

#include "../include/unit.h"

Matrix unit(const Matrix& vec){
    Matrix outvec = Matrix(3,1);
    double small, magv;
    small = 0.000001;
    magv = Matrix::norm(vec);

    if ( magv > small ) {
        for (int i = 1; i <= 3; i++) {
            outvec(i) = vec(i) / magv;
        }
    } else {
        for (int i = 1; i <= 3; i++) {
            outvec(i) = 0.0;
        }
    }
    return outvec;
}