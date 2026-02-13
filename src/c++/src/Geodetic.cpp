/// @file Geodetic.cpp
/// @brief Source file for Geodetic function
/// @date 2025-04-16
/// @author Martín Hernández Tonzán

#include <stdexcept>
#include "../include/Geodetic.h"
#include "../include/SAT_Const.h"

void Geodetic(const Matrix& r, double &lon, double &lat, double &h){
    double R_equ, f, epsRequ, e2, X, Y, Z, rho2, ZdZ,Nh, SinPhi, dZ, N, dZ_new;
    R_equ = R_Earth;
    f     = f_Earth;

    epsRequ = std::numeric_limits<double>::epsilon()*R_equ;        // Convergence criterion
    e2      = f*(2.0-f);        // Square of eccentricity

    X = r(1);                   // Cartesian coordinates
    Y = r(2);
    Z = r(3);
    rho2 = X*X + Y*Y;           // Square of distance from z-axis

                                                            // Check validity of input data
    if (Matrix::norm(r)==0.0) {
        lon = 0.0;
        lat = 0.0;
        h = -R_Earth;
        throw std::out_of_range("invalid input in Geodetic constructor\n");
    }
    // Iteration
    dZ = e2*Z;

    while(true){
        ZdZ    =  Z + dZ;
        Nh     =  sqrt ( rho2 + ZdZ*ZdZ );
        SinPhi =  ZdZ / Nh;                    // Sine of geodetic latitude
        N      =  R_equ / sqrt(1.0-e2*SinPhi*SinPhi);
        dZ_new =  N*e2*SinPhi;
        if ( fabs(dZ-dZ_new) < epsRequ )
            break;
        dZ = dZ_new;
    }

    // Longitude, latitude, altitude
    lon = atan2 ( Y, X );
    lat = atan2 ( ZdZ, sqrt(rho2) );
    h   = Nh - N;
}