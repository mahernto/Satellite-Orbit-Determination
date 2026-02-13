/// @file gast.cpp
/// @brief Source file for gast function
/// @date 2025-04-16
/// @author Martín Hernández Tonzán

#include "../include/gast.h"
#include "../include/EqnEquinox.h"
#include "../include/gmst.h"

double gast(double Mjd_UT1){
    return fmod ( gmst(Mjd_UT1) + EqnEquinox(Mjd_UT1), 2*M_PI );
}