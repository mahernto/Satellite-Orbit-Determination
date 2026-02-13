/// @file elements.h
/// @brief Header file for elements function
/// @date 2025-05-07
/// @author Martín Hernández Tonzán

#ifndef C___ELEMENTS_H
#define C___ELEMENTS_H

#include "Matrix.h"


/// @brief Computes osculating Keplerian elements from satellite state vector for elliptic orbits.
/// @param[in] y State vector (6x1): [x, y, z, vx, vy, vz]
/// @param[out] p Semilatus rectum [m]
/// @param[out] a Semimajor axis [m]
/// @param[out] e Eccentricity
/// @param[out] i Inclination [rad]
/// @param[out] Omega Longitude of the ascending node [rad]
/// @param[out] omega Argument of pericenter [rad]
/// @param[out] M Mean anomaly [rad]
void elements (const double y[6], double& p, double& a, double& e, double& i, double& Omega, double& omega, double& M);

#endif //C___ELEMENTS_H
