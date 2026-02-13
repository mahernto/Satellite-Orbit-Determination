/// @file MeanObliquity.cpp
/// @brief Source file for MeanObliquity function
/// @date 2025-04-14
/// @author Martín Hernández Tonzán

#include "../include/MeanObliquity.h"

double MeanObliquity(double Mjd_TT){
    double T = (Mjd_TT-MJD_J2000)/36525;

    return Rad *( 84381.448/3600-(46.8150+(0.00059-0.001813*T)*T)*T/3600 );
}