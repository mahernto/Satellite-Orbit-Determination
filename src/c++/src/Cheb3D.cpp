/// @file Cheb3D.cpp
/// @brief Source file for Cheb3D function
/// @date 2025-04-14
/// @author Martín Hernández Tonzán

#include "../include/Cheb3D.h"

Matrix Cheb3D(double t, int N, double Ta, double Tb, const Matrix& Cx, const Matrix& Cy, const Matrix& Cz){

    // Check validity
    if ( (t<Ta) || (Tb<t) )
        throw std::runtime_error("ERROR: Time out of range in Cheb3D::Value");

    // Clenshaw algorithm
    double tau = (2*t-Ta-Tb)/(Tb-Ta);

    Matrix f1(3,1), f2(3,1), old_f1(3,1), aux(3,1);

    for (int i = N; i >= 2; i--) {
        old_f1 = f1;
        aux(1) = Cx(i);
        aux(2) = Cy(i);
        aux(3) = Cz(i);
        f1 = 2.0 * tau * f1 - f2 + aux;
        f2 = old_f1;
    }

    aux(1) = Cx(1);
    aux(2) = Cy(1);
    aux(3) = Cz(1);

    return tau * f1 - f2 + aux;
}