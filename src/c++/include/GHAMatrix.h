/// @file GHAMatrix.h
/// @brief Header file for GHAMatrix function
/// @date 2025-04-16
/// @author Martín Hernández Tonzán

#ifndef C___GHAMATRIX_H
#define C___GHAMATRIX_H

#include "Matrix.h"

/// @brief Transformation from true equator and equinox to Earth equator and Greenwich meridian system
/// @param[in] Mjd_UT1 Modified Julian Date UT1
/// @return GHAmat Greenwich Hour Angle matrix
Matrix GHAMatrix(double Mjd_UT1);

#endif //C___GHAMATRIX_H
