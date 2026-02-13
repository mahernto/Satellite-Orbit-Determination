/// @file Accel.cpp
/// @brief Source file for Accel function
/// @date 2025-04-18
/// @author Martín Hernández Tonzán

#include "../include/Accel.h"
#include "../include/timediff.h"
#include "../include/SAT_Const.h"
#include "../include/JPL_Eph_DE430.h"
#include "../include/Mjday_TDB.h"
#include "../EKF_Global.h"
#include "../include/IERS.h"
#include "../include/PrecMatrix.h"
#include "../include/NutMatrix.h"
#include "../include/PoleMatrix.h"
#include "../include/GHAMatrix.h"
#include "../include/AccelHarmonic.h"
#include "../include/AccelPointMass.h"

Matrix Accel(double x, const Matrix& Y){

    double x_pole,y_pole,UT1_UTC,LOD,dpsi,deps,dx_pole,dy_pole,TAI_UTC;
    double UT1_TAI,UTC_GPS,UT1_GPS,TT_UTC,GPS_UTC;
    Matrix r_Mercury(3, 1), r_Venus(3, 1), r_Earth(3, 1), r_Mars(3, 1),
            r_Jupiter(3, 1), r_Saturn(3, 1), r_Uranus(3, 1), r_Neptune(3, 1),
            r_Pluto(3, 1), r_Moon(3, 1), r_Sun(3, 1);

    double Mjd_UT1, Mjd_TT, MJD_TDB;

    IERS(eop, auxParam.Mjd_UTC + x/86400, x_pole, y_pole, UT1_UTC, LOD, dpsi, deps, dx_pole, dy_pole, TAI_UTC,'l');
    timediff(UT1_UTC,TAI_UTC, UT1_TAI,UTC_GPS,UT1_GPS,TT_UTC,GPS_UTC);
    Mjd_UT1 = auxParam.Mjd_UTC + x/86400 + UT1_UTC/86400;
    Mjd_TT = auxParam.Mjd_UTC + x/86400 + TT_UTC/86400;

    Matrix P = PrecMatrix(MJD_J2000,Mjd_TT);
    Matrix N = NutMatrix(Mjd_TT);
    Matrix T = N * P;
    Matrix E = PoleMatrix(x_pole,y_pole) * GHAMatrix(Mjd_UT1) * T;

    MJD_TDB = Mjday_TDB(Mjd_TT);
    JPL_Eph_DE430(MJD_TDB, r_Mercury,r_Venus,r_Earth,r_Mars,r_Jupiter,r_Saturn,r_Uranus, r_Neptune,r_Pluto,r_Moon,r_Sun);

    // Acceleration due to harmonic gravity field
    Matrix Yaux = Y.subvector(1,3).transpose();

    Matrix a = AccelHarmonic(Yaux, E, auxParam.n, auxParam.m);

    // Luni-solar perturbations
    if (auxParam.sun){
        a = a + AccelPointMass(Yaux,r_Sun,GM_Sun);
    }

    if (auxParam.moon)
        a = a + AccelPointMass(Yaux,r_Moon,GM_Moon);

    // Planetary perturbations
    if (auxParam.planets){
        a = a + AccelPointMass(Yaux,r_Mercury,GM_Mercury);
        a = a + AccelPointMass(Yaux,r_Venus,GM_Venus);
        a = a + AccelPointMass(Yaux,r_Mars,GM_Mars);
        a = a + AccelPointMass(Yaux,r_Jupiter,GM_Jupiter);
        a = a + AccelPointMass(Yaux,r_Saturn,GM_Saturn);
        a = a + AccelPointMass(Yaux,r_Uranus,GM_Uranus);
        a = a + AccelPointMass(Yaux,r_Neptune,GM_Neptune);
        a = a + AccelPointMass(Yaux,r_Pluto,GM_Pluto);
    }
    return Y.subvector(4,6).concatRow(a.transpose()).transpose(); // return column vector
}