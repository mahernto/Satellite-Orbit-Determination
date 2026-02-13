/// @file PrecMatrix.h
/// @brief Header file for PrecMatrix function
/// @date 2025-04-16
/// @author Martín Hernández Tonzán

#ifndef C___PRECMATRIX_H
#define C___PRECMATRIX_H

#include "Matrix.h"

/// @brief Calculates precession transformation of equatorial coordinates
/// @param[in] Mjd_1 Epoch given (Modified Julian Date TT)
/// @param[in] MjD_2 Epoch to precess to (Modified Julian Date TT)
/// @return PrecMat Precession transformation matrix
Matrix PrecMatrix(double Mjd_1, double Mjd_2);

#endif //C___PRECMATRIX_H
