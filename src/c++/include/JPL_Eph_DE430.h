/// @file JPL_Eph_DE430.h
/// @brief header file for JPL_Eph_DE430 function
/// @date 2025-04-16
/// @author Martín Hernández Tonzán

#ifndef C___JPL_EPH_DE430_H
#define C___JPL_EPH_DE430_H

#include "Matrix.h"

void JPL_Eph_DE430(double Mjd_TDB,
                   Matrix& r_Mercury, Matrix& r_Venus, Matrix& r_Earth, Matrix& r_Mars,
                   Matrix& r_Jupiter, Matrix& r_Saturn, Matrix& r_Uranus, Matrix& r_Neptune,
                   Matrix& r_Pluto, Matrix& r_Moon, Matrix& r_Sun);

#endif //C___JPL_EPH_DE430_H
