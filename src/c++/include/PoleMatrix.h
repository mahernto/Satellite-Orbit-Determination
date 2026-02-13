/// @file PoleMatrix.h
/// @brief Header file for PoleMatrix function
/// @date 2025-04-16
/// @author Martín Hernández Tonzán

#ifndef C___POLEMATRIX_H
#define C___POLEMATRIX_H

#include "Matrix.h"

/// @brief Transformation from pseudo Earth-fixed to Earth-fixed coordinates for a given date
/// @param[in] xp x pole coordinate
/// @param[in] yp y pole coordinate
/// @return PoleMat Pole matrix
Matrix PoleMatrix(double xp, double yp);

#endif //C___POLEMATRIX_H
