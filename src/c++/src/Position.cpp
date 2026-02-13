/// @file Position.cpp
/// @brief Source file for Position function
/// @date 2025-04-10
/// @author Martín Hernández Tonzán

#include "../include/Position.h"
#include <cmath>


Matrix Position(double lon, double lat, double h){
    Matrix r = Matrix(3,1);
    double R_equ, f, e2, CosLat, SinLat, N;

    R_equ = R_Earth;
    f = f_Earth;

    e2 = f*(2.0-f);   // Square of eccentricity
    CosLat = cos(lat);    // (Co)sine of geodetic latitude
    SinLat = sin(lat);

    // Position vector
    N = R_equ / sqrt(1.0-e2*SinLat*SinLat);

    r(1) =  (         N+h)*CosLat*cos(lon);
    r(2) =  (         N+h)*CosLat*sin(lon);
    r(3) =  ((1.0-e2)*N+h)*SinLat;
    return r;
}