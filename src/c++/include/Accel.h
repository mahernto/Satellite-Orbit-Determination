/// @file Accel.h
/// @brief Header file for Accel function
/// @date 2025-04-18
/// @author Martín Hernández Tonzán

#ifndef C___ACCEL_H
#define C___ACCEL_H

#include "Matrix.h"

/// @ brief Computes the acceleration of an Earth orbiting satellite due to
///    - the Earth's harmonic gravity field,
///    - the gravitational perturbations of the Sun and Moon
///    - the solar radiation pressure and
///    - the atmospheric drag
/// @param[in] Mjd_TT Terrestrial Time (Modified Julian Date)
/// @param[in] Y Satellite state vector in the ICRF/EME2000 system
/// @return Matrix (6x1) dY Acceleration (a=d^2r/dt^2) in the ICRF/EME2000 system
Matrix Accel(double Mjd_TT, const Matrix& Y);

#endif //C___ACCEL_H
