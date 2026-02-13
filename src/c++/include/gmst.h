/// @file gmst.h
/// @brief Header file for gmst function
/// @date 2025-04-16
/// @author Martín Hernández Tonzán

#ifndef C___GMST_H
#define C___GMST_H

#include "SAT_Const.h"

/// @brief Computes Greenwich Mean Sidereal Time
/// @param[in] Mjd_UT1 Modified Julian Date UT1
/// @return GMST in [rad]
double gmst(double Mjd_UT1);

#endif //C___GMST_H
