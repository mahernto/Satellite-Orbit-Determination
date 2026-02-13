/// @file AccelPointMass.h
/// @brief Header file for AccelPointMass function
/// @date 2025-04-14
/// @author Martín Hernández Tonzán

#ifndef C___ACCELPOINTMASS_H
#define C___ACCELPOINTMASS_H

#include "../include/SAT_Const.h"
#include "Matrix.h"
#include <cmath>

/// @ brief Computes the pertubational acceleration due to a point mass
/// @param[in] r Satellite position vector
/// @param[in] s Point mass position vector
/// @param[in] GM Gravitational coefficient of point mass
/// @return Matrix (3x1) representing Acceleration vector (a=d^2r/dt^2)
Matrix AccelPointMass(Matrix r, const Matrix& s, double GM);

#endif //C___ACCELPOINTMASS_H
