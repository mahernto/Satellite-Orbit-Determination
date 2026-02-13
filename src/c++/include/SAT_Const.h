/// @file SAT_Const.h
/// @brief Header file for satellite constants
/// @date 2025-04-09
/// @author Martín Hernández Tonzán

#ifndef _SATCONST_
#define _SATCONST_

#include <cmath>

// Mathematical constants
const double pi2 = 2 * M_PI;              // 2pi
const double Rad = M_PI / 180.0;          // Radians per degree
const double Deg = 180.0 / M_PI;          // Degrees per radian
const double Arcs = 3600.0 * 180.0 / M_PI; // Arcseconds per radian

// General constants
const double MJD_J2000 = 51544.5;        // Modified Julian Date of J2000
const double T_B1950 = -0.500002108;     // Epoch B1950
const double c_light = 299792458.0;      // Speed of light [m/s]; DE430
const double AU = 149597870700.0;        // Astronomical unit [m]; DE430

// Physical parameters of the Earth, Sun, and Moon
const double R_Earth = 6378.1363e3;      // Earth's radius [m]; DE430
const double f_Earth = 1.0 / 298.257223563; // Flattening; WGS-84
const double R_Sun = 696000e3;           // Sun's radius [m]; DE430
const double R_Moon = 1738e3;            // Moon's radius [m]; DE430

// Earth rotation (derivative of GMST at J2000; differs from inertial period by precession)
const double omega_Earth = 15.04106717866910 / 3600.0 * Rad; // [rad/s]; WGS-84

// Gravitational coefficients
const double GM_Earth = 398600.435436e9;      // [m^3/s^2]; DE430
const double GM_Sun = 132712440041.939400e9;  // [m^3/s^2]; DE430
const double GM_Moon = GM_Earth / 81.30056907419062; // [m^3/s^2]; DE430
const double GM_Mercury = 22031.780000e9;     // [m^3/s^2]; DE430
const double GM_Venus = 324858.592000e9;      // [m^3/s^2]; DE430
const double GM_Mars = 42828.375214e9;        // [m^3/s^2]; DE430
const double GM_Jupiter = 126712764.800000e9; // [m^3/s^2]; DE430
const double GM_Saturn = 37940585.200000e9;   // [m^3/s^2]; DE430
const double GM_Uranus = 5794548.600000e9;    // [m^3/s^2]; DE430
const double GM_Neptune = 6836527.100580e9;   // [m^3/s^2]; DE430
const double GM_Pluto = 977.0000000000009e9;  // [m^3/s^2]; DE430

// Solar radiation pressure at 1 AU
const double P_Sol = 1367 / c_light;  // [N/m^2] (~1367 W/m^2); IERS 96

#endif

