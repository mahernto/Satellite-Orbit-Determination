/// @file NutAngles.h
/// @brief Header file for NutAngles function
/// @date 2025-04-10
/// @author Martín Hernández Tonzán

#ifndef C___NUTANGLES_H
#define C___NUTANGLES_H

#include "../include/SAT_Const.h"
#include <utility>

/// @ brief Nutation in longitude and obliquity
/// @param[in] Mjd_TT Modified julian date (Terrestrial Time)
/// @return pair containing Nutation angles dpsi & deps
void NutAngles(double Mjd_TT, double &dpsi, double &deps);

#endif //C___NUTANGLES_H
