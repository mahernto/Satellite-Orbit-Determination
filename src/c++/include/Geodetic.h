/// @file Geodetic.h
/// @brief Header file for Geodetic function
/// @date 2025-04-16
/// @author Martín Hernández Tonzán

#ifndef C___GEODETIC_H
#define C___GEODETIC_H

#include "Matrix.h"
#include <limits>

/// @brief Computes geodetic coordinates (Longitude [rad], latitude [rad], altitude [m])
/// @param[in] r position vector [m]
/// @param[out] lon longitude [rad]
/// @param[out] lat latitude [rad]
/// @param[out] h altitude [m]
void Geodetic(const Matrix& r, double &lon, double &lat, double &h);

#endif //C___GEODETIC_H
