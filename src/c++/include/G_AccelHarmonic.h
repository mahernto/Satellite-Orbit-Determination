/// @file G_AccelHarmonic.h
/// @brief Header file for G_AccelHarmonic function
/// @date 2025-04-16
/// @author Martín Hernández Tonzán


#ifndef C___G_ACCELHARMONIC_H
#define C___G_ACCELHARMONIC_H

#include "Matrix.h"

/// @ brief Computes the gradient of the Earth's harmonic gravity field
/// @param[in] r Satellite position vector in the true-of-date system
/// @param[in] E Transformation matrix to body-fixed system
/// @param[in] n_max Gravity model degree
/// @param[in] m_max Gravity model order
/// @return Matrix representing Gradient (G=da/dr) in the true-of-date system
Matrix G_AccelHarmonic(Matrix r, const Matrix& U, int n_max, int m_max);

#endif //C___G_ACCELHARMONIC_H
