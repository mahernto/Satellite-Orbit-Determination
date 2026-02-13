/// @file G_AccelHarmonic.cpp
/// @brief Source file for G_AccelHarmonic function
/// @date 2025-04-16
/// @author Martín Hernández Tonzán

#include "../include/G_AccelHarmonic.h"
#include "../include/AccelHarmonic.h"
#include "../EKF_Global.h"

Matrix G_AccelHarmonic(Matrix r, const Matrix& U, int n_max, int m_max){
    double d = 1.0;   // Position increment [m]

    Matrix G = Matrix(3,3);
    Matrix dr = Matrix(3,1);

    // Gradient
    for (int i=1; i <= 3; i++){
        // Set offset in i-th component of the position vector
        dr = Matrix(3,1);
        dr(i) = d;
        // Acceleration difference
        Matrix da = AccelHarmonic ( r+dr/2,U, n_max, m_max) - AccelHarmonic ( r-dr/2,U, n_max, m_max);
        // Derivative with respect to i-th axis
        G(1,i) = (da/d)(1);
        G(2,i) = (da/d)(2);
        G(3,i) = (da/d)(3);
    }
    return G;
}