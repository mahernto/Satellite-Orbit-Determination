/// @file NutMatrix.h
/// @brief Header file for NutMatrix function
/// @date 2025-04-16
/// @author Martín Hernández Tonzán

#ifndef C___NUTMATRIX_H
#define C___NUTMATRIX_H

#include "Matrix.h"

/// @brief Transformation from mean to true equator and equinox
/// @param[in] Mjd_TT Modified Julian Date (Terrestrial Time)
/// @return NutMat Nutation matrix
Matrix NutMatrix(double Mjd_TT);

#endif //C___NUTMATRIX_H
