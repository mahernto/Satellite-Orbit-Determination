//
// Created by martin on 23/04/2025.
//

#ifndef C___DEINTEG_H
#define C___DEINTEG_H

// y = DEInteg(func,t,tout,relerr,abserr,n_eqn,y)

#include "Matrix.h"
#include <limits>
#include <cmath>

/// @ brief Integrates a system of differential equations using a numerical method
/// @param[in] func Pointer to the function that evaluates the derivatives of the system function must have the following form: Matrix func(double t, const Matrix& y)
/// @param[in] t Initial time
/// @param[int] tout Target time to integrate to
/// @param[in] relerr Relative error tolerance for adaptive step size control
/// @param[in] abserr Absolute error tolerance for adaptive step size control
/// @param[in] n_eqn Number of equations (dimension of the system)
/// @param[in,out] y On input, the initial state vector (size n_eqn x 1); On output, the state vector at `tout`
void DEInteg(Matrix (*func)(double, const Matrix &), double t, double tout, double relerr, double abserr, int n_eqn, Matrix& y);

#endif //C___DEINTEG_H
