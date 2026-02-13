/// @file anglesg.h
/// @brief Header file for anglesg function
/// @date 2025-05-07
/// @author Martín Hernández Tonzán

#ifndef C___ANGLESG_H
#define C___ANGLESG_H

#include "Matrix.h"

/// @brief Orbit determination using three optical angle-only observations.
/// @param[in] az1   Azimuth at observation time t1 [rad]
/// @param[in] az2   Azimuth at observation time t2 [rad]
/// @param[in] az3   Azimuth at observation time t3 [rad]
/// @param[in] el1   Elevation at observation time t1 [rad]
/// @param[in] el2   Elevation at observation time t2 [rad]
/// @param[in] el3   Elevation at observation time t3 [rad]
/// @param[in] Mjd1  Modified Julian Date at time t1
/// @param[in] Mjd2  Modified Julian Date at time t2
/// @param[in] Mjd3  Modified Julian Date at time t3
/// @param[in] Rs1   IJK position vector of observation site at t1 [m]
/// @param[in] Rs2   IJK position vector of observation site at t2 [m]
/// @param[in] Rs3   IJK position vector of observation site at t3 [m]
/// @param[out] r2   IJK position vector of the satellite at t2 [m]
/// @param[out] v2   IJK velocity vector of the satellite at t2 [m/s]
void anglesg(double az1, double az2, double az3, double el1,double el2, double el3, double Mjd1, double Mjd2, double Mjd3, Matrix Rs1, Matrix Rs2, Matrix Rs3, Matrix &r2, Matrix &v2);

#endif //C___ANGLESG_H
