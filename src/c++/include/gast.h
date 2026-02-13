/// @file gast.h
/// @brief Header file for gast function
/// @date 2025-04-16
/// @author Martín Hernández Tonzán

#ifndef C___GAST_H
#define C___GAST_H

/// @brief Computes Greenwich Apparent Sidereal Time
/// @param[in] Mjd_UT1 Modified Julian Date UT1
/// @return GAST in [rad]
double gast(double Mjd_UT1);

#endif //C___GAST_H
