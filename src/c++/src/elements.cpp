/// @file elements.cpp
/// @brief Source file for elements function
/// @date 2025-05-07
/// @author Martín Hernández Tonzán

#include "../include/elements.h"
#include "../include/SAT_Const.h"
#include "vector.h"

void elements (const double y[6], double& p, double& a, double& e, double& i, double& Omega, double& omega, double& M){

    double r[3] = { y[0], y[1], y[2] };                // Position
    double v[3] = { y[3], y[4], y[5] };                // Velocity

    double h[3];
    cross(r, v, h);                                    // Areal velocity
    double magh = norm(h);
    p = magh * magh / GM_Earth;

    double H = magh;

    Omega = std::atan2(h[0], -h[1]);                         // Long. ascend. node
    Omega = std::fmod(Omega + pi2, pi2);

    i = std::atan2(std::sqrt(h[0]*h[0] + h[1]*h[1]), h[2]); // Inclination
    double u = std::atan2(r[2]*H, -r[0]*h[1] + r[1]*h[0]);     // Arg. of latitude

    double R = norm(r);                                           // Distance

    a = 1.0 / (2.0 / R - dot(v, v) / GM_Earth);                // Semi-major axis

    double eCosE = 1.0 - R / a;                                      // e*cos(E)
    double eSinE = dot(r, v) / std::sqrt(GM_Earth * a);     // e*sin(E)

    double e2 = eCosE*eCosE + eSinE*eSinE;
    e = std::sqrt(e2);                                            // Eccentricity
    double E = std::atan2(eSinE, eCosE);                       // Eccentric anomaly

    M = std::fmod(E - eSinE + pi2, pi2);                       // Mean anomaly

    double nu = std::atan2(std::sqrt(1.0 - e2) * eSinE, eCosE - e2); // True anomaly

    omega = std::fmod(u - nu + pi2, pi2);                       // Arg. of perihelion
}
