/// @file IERS.h
/// @brief Header file for IERS function
/// @date 2025-04-16
/// @author Martín Hernández Tonzán

#ifndef C___IERS_H
#define C___IERS_H

#include "Matrix.h"
#include "SAT_Const.h"

/// @brief Interpolates Earth orientation parameters from IERS data
/// @param[in] eop Matrix containing Earth orientation parameters
/// @param[in] Mjd_UTC Modified Julian Date in UTC
/// @param[out] x_pole Polar motion in the x direction
/// @param[out] y_pole Polar motion in the y direction
/// @param[out] UT1_UTC Difference UT1 - UTC
/// @param[out] LOD Length of day variation
/// @param[out] dpsi Nutation correction in longitude
/// @param[out] deps Nutation correction in obliquity
/// @param[out] dx_pole Celestial pole offset in x
/// @param[out] dy_pole Celestial pole offset in y
/// @param[out] TAI_UTC Difference TAI - UTC
/// @param[in] interp Interpolation method: 'l' for linear, 'n' for nearest-neighbor
void IERS(const Matrix& eop, double Mjd_UTC, double& x_pole, double& y_pole, double& UT1_UTC, double& LOD, double& dpsi, double& deps, double& dx_pole, double& dy_pole, double& TAI_UTC, char interp='n');

#endif //C___IERS_H
