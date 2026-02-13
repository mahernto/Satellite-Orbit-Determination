/// @file EqnEquinox.cpp
/// @brief Source file for EqnEquinox function
/// @date 2025-04-16
/// @author Martín Hernández Tonzán

#include "../include/EqnEquinox.h"
#include "../include/NutAngles.h"
#include "../include/MeanObliquity.h"

double EqnEquinox(double Mjd_TT){
    double dpsi, deps;
    NutAngles(Mjd_TT, dpsi, deps);
    return dpsi * cos(MeanObliquity(Mjd_TT));
}