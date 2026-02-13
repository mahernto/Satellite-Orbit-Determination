/// @file Legendre.h
/// @brief Header file for Legendre function
/// @date 2025-04-10
/// @author Martín Hernández Tonzán

#ifndef C___LEGENDRE_H
#define C___LEGENDRE_H

#include <tuple>
#include "Matrix.h"

/// @brief Solves Legendre polynomials
/// @param[in] n maximum degree
/// @param[in] m maximum order
/// @param[in] fi angle in radians
/// @param[out] pnm Matrix which stores computed Legendre polynomials
/// @param[out] dpnm Matrix which stores the derivatives with respect to fi
void Legendre(double n, double m, double fi, Matrix &pnm, Matrix &dpnm);

#endif //C___LEGENDRE_H
