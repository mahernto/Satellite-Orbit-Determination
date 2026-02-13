/// @file LTC.h
/// @brief Header file for LTC function
/// @date 2025-04-16
/// @author Martín Hernández Tonzán

#ifndef C___LTC_H
#define C___LTC_H

#include "Matrix.h"

/// @brief Transformation from Greenwich meridian system to local tangent coordinates
/// @param[in] lon Geodetic East longitude [rad]
/// @param[in] lat Geodetic latitude [rad]
/// @return M Rotation matrix from the Earth equator and Greenwich meridian to the local tangent (East-North-Zenith) coordinate system
Matrix LTC(double lon, double lat);

#endif //C___LTC_H
