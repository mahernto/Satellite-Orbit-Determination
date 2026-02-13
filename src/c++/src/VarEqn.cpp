/// @file VarEqn.cpp
/// @brief Source file for VarEqn function
/// @date 2025-04-16
/// @author Martín Hernández Tonzán

#include "../include/VarEqn.h"
#include "../EKF_Global.h"
#include "../include/IERS.h"
#include "../include/timediff.h"
#include "../include/PrecMatrix.h"
#include "../include/NutMatrix.h"
#include "../include/PoleMatrix.h"
#include "../include/GHAMatrix.h"
#include "../include/AccelHarmonic.h"
#include "../include/G_AccelHarmonic.h"

Matrix VarEqn(double x, const Matrix& yPhi) {
    double x_pole, y_pole, UT1_UTC, LOD, dpsi, deps, dx_pole, dy_pole, TAI_UTC;
    double UT1_TAI, UTC_GPS, UT1_GPS, TT_UTC, GPS_UTC;

    IERS(eop, auxParam.Mjd_UTC, x_pole, y_pole, UT1_UTC, LOD, dpsi, deps, dx_pole, dy_pole, TAI_UTC, 'l');

    timediff(UT1_UTC, TAI_UTC, UT1_TAI, UTC_GPS, UT1_GPS, TT_UTC, GPS_UTC);
    double Mjd_UT1 = auxParam.Mjd_TT + (UT1_UTC - TT_UTC) / 86400;

    // Transformation matrix
    Matrix P = PrecMatrix(MJD_J2000, auxParam.Mjd_TT + x / 86400);
    Matrix N = NutMatrix(auxParam.Mjd_TT + x / 86400);
    Matrix T = N * P;
    Matrix E = PoleMatrix(x_pole, y_pole) * GHAMatrix(Mjd_UT1) * T;

    // State vector components
    Matrix r = yPhi.subvector(1, 3);
    Matrix v = yPhi.subvector(4, 6);
    Matrix Phi(6, 6);

    // State transition matrix
    for (int j = 1; j <= 6; j++) {
        for (int i = 1; i <= 6; ++i) {
            Phi(i, j) = yPhi(6 * j + i);
        }
    }

    // Acceleration and gradient
    Matrix a = AccelHarmonic(r.transpose(), E, auxParam.n, auxParam.m);
    Matrix G = G_AccelHarmonic(r.transpose(), E, auxParam.n, auxParam.m);

    // Time derivative of state transition matrix
    Matrix yPhip(42, 1);
    Matrix dfdy(6, 6);

    for (int i = 1; i <= 3; i++) {
        for (int j = 1; j <= 3; j++) {
            dfdy(i, j) = 0.0;                 // dv/dr(i,j)
            dfdy(i + 3, j) = G(i, j);            // da/dr(i,j)
            if (i == j){
                dfdy(i, j + 3) = 1;
            }
            else{
                dfdy(i, j + 3) = 0;            // dv/dv(i,j)
            }
            dfdy(i + 3, j + 3) = 0.0;         // da/dv(i,j)
        }
    }

    Matrix Phip = dfdy * Phi;

    // Derivative of combined state vector and state transition matrix
    for (int i = 1; i <= 3; i++) {
        yPhip(i) = v(i);                 // dr/dt(i)
        yPhip(i + 3) = a(i);                 // dv/dt(i)
    }
    for (int i = 1; i <= 6; i++) {
        for (int j = 1; j <= 6; j++) {
            yPhip(6 * j + i) = Phip(i, j);     // dPhi/dt(i,j)
        }
    }

    return yPhip;
}