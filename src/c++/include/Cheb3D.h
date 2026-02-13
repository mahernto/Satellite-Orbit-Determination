/// @file Cheb3D.h
/// @brief Header file for Cheb3D function
/// @date 2025-04-14
/// @author Martín Hernández Tonzán

#ifndef C___CHEB3D_H
#define C___CHEB3D_H

#include <iostream>
#include "Matrix.h"

/// @ brief Chebyshev approximation of 3-dimensional vectors
/// @param[in] t The time at which the Chebyshev polynomial should be evaluated.
/// @param[in] N Number of coefficients
/// @param[in] Ta Begin interval
/// @param[in] Tb End interval
/// @param[in] Cx Coefficients of Chebyshev polyomial (x-coordinate)
/// @param[in] Cy Coefficients of Chebyshev polyomial (y-coordinate)
/// @param[in] Cz Coefficients of Chebyshev polyomial (z-coordinate)
/// @return Matrix (3x1) representing evaluated vector at time t
Matrix Cheb3D(double t, int N, double Ta, double Tb, const Matrix& Cx, const Matrix& Cy, const Matrix& Cz);

#endif //C___CHEB3D_H
