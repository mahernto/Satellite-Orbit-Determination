/// @file Position.h
/// @brief Header file for Position function
/// @date 2025-04-10
/// @author Martín Hernández Tonzán

#ifndef C___POSITION_H
#define C___POSITION_H

#include "./SAT_Const.h"
#include "Matrix.h"

/// @brief Position vector from geodetic coordinates (Longitude [rad], latitude [rad], altitude [m])
/// @param[in] lon longitude
/// @param[in] lat latitude
/// @param[in] h altitude
/// @return Matrix representing result vector
Matrix Position(double lon, double lat, double h);


#endif //C___POSITION_H
