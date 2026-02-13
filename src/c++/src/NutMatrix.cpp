/// @file NutMatrix.cpp
/// @brief Source file for NutMatrix function
/// @date 2025-04-16
/// @author Martín Hernández Tonzán

#include "../include/NutMatrix.h"
#include "../include/MeanObliquity.h"
#include "../include/NutAngles.h"
#include "../include/R_z.h"
#include "../include/R_x.h"

Matrix NutMatrix(double Mjd_TT){
    // Mean obliquity of the ecliptic
    double eps = MeanObliquity (Mjd_TT);

    double dpsi, deps;

    // Nutation in longitude and obliquity
    NutAngles (Mjd_TT, dpsi, deps);

    // Transformation from mean to true equator and equinox
    return R_x(-eps-deps)*R_z(-dpsi)*R_x(+eps);
}
