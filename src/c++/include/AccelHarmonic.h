/// @file AccelHarmonic.h
/// @brief Header file for AccelHarmonic function
/// @date 2025-04-15
/// @author Martín Hernández Tonzán

#ifndef C___ACCELHARMONIC_H
#define C___ACCELHARMONIC_H

#include "Matrix.h"
#include "SAT_Const.h"

/// @ brief Computes the acceleration due to the harmonic gravity field of the central body
/// @param[in] r Satellite position vector in the inertial system
/// @param[in] E Transformation matrix to vody-fixed system
/// @param[in] n_max Maximum degree
/// @param[in] m_max Maximum order (m_max<=n_max; m_max=0 for zonals, only)
/// @return Matrix (3x1) representing Acceleration vector (a=d^2r/dt^2)
Matrix AccelHarmonic(const Matrix& r, Matrix E, int n_max, int m_max);


#endif //C___ACCELHARMONIC_H
