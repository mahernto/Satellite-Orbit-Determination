/// @file MeanObliquity.h
/// @brief Header file for MeanObliquity function
/// @date 2025-04-14
/// @author Martín Hernández Tonzán


#ifndef C___MEANOBLIQUITY_H
#define C___MEANOBLIQUITY_H
#include "../include/SAT_Const.h"

/// @ brief Computes the mean obliquity of the ecliptic
/// @param[in] Mjd_TT Modified julian date (Terrestrial Time)
/// @return MOblq Mean obliquity of the ecliptic [rad]
double MeanObliquity(double Mjd_TT);

#endif //C___MEANOBLIQUITY_H
