/// @file EKF_Test.cpp
/// @brief Source file for unit testing
/// @date 2025-04-09
/// @version 1.0

#include <cstdio>
#include <iostream>
#include <cmath>
#include <iomanip>
#include <cstring>

#include "include/SAT_Const.h"
#include "./include/Matrix.h"
#include "EKF_Global.h"

#include "include/Mjday.h"
#include "include/R_x.h"
#include "include/R_y.h"
#include "include/R_z.h"
#include "include/Legendre.h"
#include "include/Frac.h"
#include "include/Position.h"
#include "include/AccelPointMass.h"
#include "include/Mjday_TDB.h"
#include "include/NutAngles.h"
#include "include/MeanObliquity.h"
#include "include/Cheb3D.h"
#include "include/timediff.h"
#include "include/AzElPa.h"
#include "include/unit.h"
#include "include/MeasUpdate.h"
#include "include/TimeUpdate.h"
#include "include/sign_.h"
#include "include/EqnEquinox.h"
#include "include/gmst.h"
#include "include/AccelHarmonic.h"
#include "include/gast.h"
#include "include/G_AccelHarmonic.h"
#include "include/LTC.h"
#include "include/NutMatrix.h"
#include "include/PoleMatrix.h"
#include "include/PrecMatrix.h"
#include "include/GHAMatrix.h"
#include "include/Geodetic.h"
#include "include/IERS.h"
#include "include/Accel.h"
#include "include/JPL_Eph_DE430.h"
#include "include/VarEqn.h"
#include "include/DEInteg.h"
#include "include/angl.h"
#include "include/gibbs.h"
#include "include/elements.h"
#include "include/hgibbs.h"
#include "include/anglesg.h"

using namespace std;

Matrix eop(13, 21413);
Matrix PC(2285, 1020);
Matrix obs(46, 4);

Matrix Cnm = Matrix(181, 181);
Matrix Snm = Matrix(181, 181);

AuxParam auxParam = AuxParam();

void loadEOP(const char * filepath){
    FILE *fp;
    fp = fopen(filepath, "r");
    if (fp == nullptr) {
        printf("Fail open eop19620101.txt file\n");
        exit(EXIT_FAILURE);
    }

    for (int i = 1; i <= 21413; i++) {
        fscanf(fp, "%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf%lf", &eop(1, i), &eop(2, i), &eop(3, i), &eop(4, i),
               &eop(5, i), &eop(6, i), &eop(7, i), &eop(8, i), &eop(9, i), &eop(10, i), &eop(11, i), &eop(12, i),
               &eop(13, i));
    }
    fclose(fp);
}

void loadDE430Coeff(const char * filepath){
    FILE *fp;
    fp = fopen(filepath, "r");
    if (fp == nullptr) {
        printf("Fail open DE430Coeff.txt file\n");
        exit(EXIT_FAILURE);
    }

    for (int i = 1; i <= 2285; i++) {
        for (int j = 1; j <= 1020; j++) {
            fscanf(fp, "%lf", &PC(i, j));
        }
    }
    fclose(fp);
}

void loadGMS03S(const char * filepath){
    FILE *fp;
    fp = fopen(filepath, "r");
    if (fp == nullptr) {
        printf("Fail open GGM03S.txt file\n");
        exit(EXIT_FAILURE);
    }
    double f, c, aux1, aux2;
    for (int n = 0; n <= 180; n++) {
        for (int m = 0; m <= n; m++) {
            fscanf(fp, "%lf%lf%lf%lf%lf%lf", &f, &c, &Cnm(n + 1, m + 1), &Snm(n + 1, m + 1), &aux1, &aux2);
        }
    }
    fclose(fp);
}

void loadGEOS3(const char * filepath){
    FILE *fp;
    //read Earth orientation parameters
    //  ----------------------------------------------------------------------------------------------------
    // |  Date    MJD      x         y       UT1-UTC      LOD       dPsi    dEpsilon     dX        dY    DAT
    // |(0h UTC)           "         "          s          s          "        "          "         "     s
    //  ----------------------------------------------------------------------------------------------------
    // read observations
    fp = fopen(filepath, "r");
    if (fp == nullptr) {
        printf("Fail open GEOS3.txt file\n");
        exit(EXIT_FAILURE);
    }
    char line[55], y[5], mo[3], d[3], h[3], mi[3], seconds[7], a[9], e[8], di[11];
    double ss, az, el, Distance;
    int Year, M, D, hh, mm;
    for (int i = 1; i <= 46; i++) {
        fgets(line, sizeof(line) + 2, fp);

        strncpy(y, &line[0], 4);
        y[4] = '\0';
        Year = atoi(y);

        strncpy(mo, &line[5], 2);
        mo[2] = '\0';
        M = atoi(mo);

        strncpy(d, &line[8], 2);
        d[2] = '\0';
        D = atoi(d);

        strncpy(h, &line[12], 2);
        h[2] = '\0';
        hh = atoi(h);

        strncpy(mi, &line[15], 2);
        mi[2] = '\0';
        mm = atoi(mi);

        strncpy(seconds, &line[18], 6);
        seconds[6] = '\0';
        ss = atof(seconds);

        strncpy(a, &line[25], 8);
        a[8] = '\0';
        az = atof(a);

        strncpy(e, &line[35], 7);
        e[7] = '\0';
        el = atof(e);

        strncpy(di, &line[44], 10);
        di[10] = '\0';
        Distance = atof(di);

        obs(i, 1) = Mjday(Year, M, D, hh, mm, ss);
        obs(i, 2) = Rad * az;
        obs(i, 3) = Rad * el;
        obs(i, 4) = 1e3 * Distance;
    }
    fclose(fp);
}

int tests_run = 0;

#define FAIL() printf("\nfailure in %s() line %d\n", __func__, __LINE__)
#define _assert(test) do { if (!(test)) { FAIL(); return 1; } } while(0)
#define _verify(test) do { int r=test(); tests_run++; if(r) return r; } while(0)

#define TOL_ 10e-14

int Mjday_01() {
    _assert(fabs(Mjday(2025, 4, 3, 15, 37, 5) - 60768.6507523148) < pow(10, -10));
    return 0;
}

int Mjday_02() {
    _assert(fabs(Mjday(2025, 4, 3, 0, 0, 0.0) - Mjday(2025, 4, 3)) < pow(10, -10));
    return 0;
}

int R_x_01() {
    double alpha = 1.0;
    Matrix sol(3, 3);
    sol = R_x(alpha);

    _assert(fabs(sol(1, 1)) - 1 < TOL_ && fabs(sol(1, 2)) < TOL_ && fabs(sol(1, 3)) < TOL_);
    _assert(fabs(sol(2, 1)) < TOL_ && fabs(sol(2, 2) - 0.54030230586814) < TOL_ &&
            fabs(sol(2, 3) - 0.841470984807897) < TOL_);
    return 0;
}

int R_y_01() {
    double alpha = 1.0;
    Matrix sol(3, 3);
    sol = R_y(alpha);

    _assert(fabs(sol(2, 2)) - 1 < TOL_ && fabs(sol(1, 2)) < TOL_ && fabs(sol(2, 1)) < TOL_);
    _assert(fabs(sol(2, 3)) < TOL_ && fabs(sol(1, 1) - 0.54030230586814) < TOL_ &&
            fabs(sol(3, 1) - 0.841470984807897) < TOL_);
    return 0;
}

int R_z_01() {
    double alpha = 1.0;
    Matrix sol(3, 3);
    sol = R_z(alpha);

    _assert(fabs(sol(3, 3)) - 1 < TOL_ && fabs(sol(1, 3)) < TOL_ && fabs(sol(3, 1)) < TOL_);
    _assert(fabs(sol(3, 2)) < TOL_ && fabs(sol(1, 1) - 0.54030230586814) < TOL_ &&
            fabs(sol(1, 2) - 0.841470984807897) < TOL_);
    return 0;
}

int Legendre_01() {
    Matrix pnm(3, 3), dpnm(3, 3);
    Legendre(2, 2, 1, pnm, dpnm);

    _assert(fabs(pnm(1, 1)) - 1 < TOL_ && fabs(pnm(1, 2)) < TOL_ && fabs(pnm(1, 3)) < TOL_);
    _assert(fabs(pnm(2, 1)) - 1.457470498782296 < TOL_ && fabs(pnm(2, 2) - 0.935831045210238) < TOL_ &&
            fabs(pnm(2, 3)) < TOL_);
    _assert(fabs(pnm(3, 1)) - 1.256916455730625 < TOL_ && fabs(pnm(3, 2) - 1.760846895422561) < TOL_ &&
            fabs(pnm(3, 3) - 0.565313394670859) < TOL_);

    _assert(fabs(dpnm(1, 1)) < TOL_ && fabs(dpnm(1, 2)) < TOL_ && fabs(dpnm(1, 3)) < TOL_);
    _assert(fabs(dpnm(2, 1)) - 0.935831045210238 < TOL_ && fabs(dpnm(2, 2) + 1.457470498782296) < TOL_ &&
            fabs(dpnm(2, 3)) < TOL_);
    _assert(fabs(dpnm(3, 1)) - 3.049876287221798 < TOL_ && fabs(dpnm(3, 2) + 1.611729767523982) < TOL_ &&
            fabs(dpnm(3, 3) + 1.760846895422562) < TOL_);

    return 0;
}

int Frac_01() {
    _assert(fabs(Frac(1.123456789123456789)) - 0.123456789123457 < TOL_);
    return 0;
}

int Position_01() {
    Matrix sol = Position(Rad * 21.5748, Rad * (-158.2706), 300.20);
    _assert(fabs(sol(1) + 5512602.017546574) < TOL_ * 10e6 &&
            fabs(sol(2) + 2179789.386247879) < TOL_ * 10e6 &&
            fabs(sol(3) + 2346716.056112798) < TOL_ * 10e6);
    return 0;
}

int AccelPointMass_01() {
    double rarray[3] = {1.0, 8.0, 2.0};
    double sarray[3] = {2.0, 7.0, 12.0};
    Matrix r = Matrix(3, 1, rarray, 3);
    Matrix s = Matrix(3, 1, sarray, 3);
    Matrix m = AccelPointMass(r, s, GM_Sun);

    _assert(fabs(m(1) * 10e-17 - 3.28347272984212) < TOL_ &&
            fabs(m(2) * 10e-18 + 4.648059387053737) < TOL_ &&
            fabs(m(3) * 10e-18 - 7.123216831237252) < TOL_);
    return 0;
}

int AccelPointMass_02() {
    double rarray[3] = {1.0, 3.5, 6.0};
    double sarray[3] = {2.0, 3.2, 12.0};
    Matrix r = Matrix(3, 1, rarray, 3);
    Matrix s = Matrix(3, 1, sarray, 3);
    Matrix m = AccelPointMass(r, s, GM_Sun);

    _assert(fabs(m(1) * 10e-18 - 4.54182966772213) < TOL_ &&
            fabs(m(2) * 10e-18 + 3.89604837976212) < TOL_ &&
            fabs(m(3) * 10e-19 - 2.72509780063328) < TOL_);
    return 0;
}

int Mjd_TDB_01() {
    double resultado = Mjday_TDB(Mjday(2025, 4, 3, 15, 37, 5));
    _assert(fabs(resultado - 60768.65075233379) < TOL_ * 10e4);
    return 0;
}

int Mjd_TDB_02() {
    double resultado = Mjday_TDB(Mjday(2025, 12, 22, 12, 5, 0.0));
    _assert(fabs(resultado - 61031.50347221790) < TOL_ * 10e4);
    return 0;
}

int NutAngles_01() {
    double dpsi, deps;
    NutAngles(61031.50347221790, dpsi, deps);
    _assert(fabs(dpsi - 0.00002510857350142765) < TOL_ &&
            fabs(deps - 0.00003920137304829823) < TOL_);
    return 0;
}

int NutAngles_02() {
    double dpsi, deps;
    NutAngles(58000.0, dpsi, deps);
    _assert(fabs(dpsi + 0.00004703913713081135) < TOL_ &&
            fabs(deps + 0.00003367682654877015) < TOL_);
    return 0;
}

int MeanObliquity_01() {
    double mo = MeanObliquity(61031.50347221790);
    _assert(fabs(mo - 0.409033852158135) < TOL_);
    return 0;
}

int MeanObliquity_02() {
    double mo = MeanObliquity(58000.0);
    _assert(fabs(mo - 0.40905268985035) < TOL_);
    return 0;
}

int Cheb3D_01() {
    const int N = 3;
    double t = 0.5;
    double Ta = 0.0;
    double Tb = 1;

    double Cxa[N] = {1.0, 0.5, 0.2};
    double Cya[N] = {0.8, 0.4, 0.1};
    double Cza[N] = {0.6, 0.3, 0.15};

    Matrix r = Cheb3D(t, N, Ta, Tb, Matrix(N, 1, Cxa, N), Matrix(N, 1, Cya, N), Matrix(N, 1, Cza, N));

    _assert(fabs(r(1) - 0.8) < TOL_ &&
            fabs(r(2) - 0.7) < TOL_ &&
            fabs(r(3) - 0.45) < TOL_);
    return 0;
}

int Cheb3D_02() {
    const int N = 5;
    double t = 100.0;
    double Ta = 0.0;
    double Tb = 400.0;

    double Cxa[N] = {1.0, 2.0, 3.0, 4.0, 5.0};
    double Cya[N] = {1.0, 0.0, -0.5, 0.75, -0.25};
    double Cza[N] = {3.0, 2.0, 1.0, 0.0, -1.0};

    Matrix r = Cheb3D(t, N, Ta, Tb, Matrix(N, 1, Cxa, N), Matrix(N, 1, Cya, N), Matrix(N, 1, Cza, N));

    _assert(fabs(r(1) - 0.0) < TOL_ &&
            fabs(r(2) - 2.125) < TOL_ &&
            fabs(r(3) - 2.0) < TOL_);
    return 0;
}

int timediff_01() {
    double UT1_UTC = 0.5;
    double TAI_UTC = 37.0;

    double UT1_TAI, UTC_GPS, UT1_GPS, TT_UTC, GPS_UTC;

    timediff(UT1_UTC, TAI_UTC, UT1_TAI, UTC_GPS, UT1_GPS, TT_UTC, GPS_UTC);

    _assert(fabs(UT1_TAI + 36.5) < TOL_ &&
            fabs(UTC_GPS + 18) < TOL_ &&
            fabs(UT1_GPS + 17.5) < TOL_ &&
            fabs(TT_UTC - 69.184) < TOL_ &&
            fabs(GPS_UTC - 18) < TOL_);
    return 0;
}

int timediff_02() {
    double UT1_UTC = 45.5;
    double TAI_UTC = 900.0;

    double UT1_TAI, UTC_GPS, UT1_GPS, TT_UTC, GPS_UTC;

    timediff(UT1_UTC, TAI_UTC, UT1_TAI, UTC_GPS, UT1_GPS, TT_UTC, GPS_UTC);

    _assert(fabs(UT1_TAI + 854.5) < TOL_ &&
            fabs(UTC_GPS + 881) < TOL_ &&
            fabs(UT1_GPS + 835.5) < TOL_ &&
            fabs(TT_UTC - 932.184) < TOL_ &&
            fabs(GPS_UTC - 881) < TOL_);
    return 0;
}

int AzElPa_01() {
    double sarray[3] = {50.0, 100.0, 500.0};

    Matrix s = Matrix(3, 1, sarray, 3);

    double Az, El;
    Matrix dAds(3, 1);
    Matrix dEds(3, 1);

    AzElPa(s, Az, El, dAds, dEds);

    _assert(fabs(Az - 0.463647609000806) < TOL_ &&
            fabs(El - 1.35080834939944) < TOL_ &&
            fabs(dAds(1) - 0.008) < TOL_ &&
            fabs(dAds(2) + 0.004) < TOL_ &&
            fabs(dAds(3) - 0) < TOL_ &&
            fabs(dEds(1) + 0.00085183541999992) < TOL_ &&
            fabs(dEds(2) + 0.00170367083999984) < TOL_ &&
            fabs(dEds(3) - 0.00042591770999996) < TOL_);

    return 0;
}

int unit_01() {
    double array[3] = {1.0, 3.0, 5.0};
    Matrix in = Matrix(3, 1, array, 3);

    Matrix out = unit(in);
    _assert(fabs(out(1) - 0.169030850945703) < TOL_ &&
            fabs(out(2) - 0.507092552837110) < TOL_ &&
            fabs(out(3) - 0.845154254728517) < TOL_);
    return 0;
}

int unit_02() {
    double array[3] = {1.0, 3.0, 500.0};
    Matrix in = Matrix(3, 1, array, 3);

    Matrix out = unit(in);
    _assert(fabs(out(1) - 0.001999960001200) < TOL_ &&
            fabs(out(2) - 0.005999880003600) < TOL_ &&
            fabs(out(3) - 0.999980000599980) < TOL_);
    return 0;
}

int MeasUpdate_01() {
    double x_data[3] = {1.0, 2.0, 3.0};
    double z_data[3] = {1.1, 2.1, 2.9};
    double g_data[3] = {1.0, 2.0, 3.0};
    double s_data[3] = {0.1, 0.1, 0.1};
    double G_data[9] = {1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0};
    double P_data[9] = {1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0};

    Matrix x = Matrix(3, 1, x_data, 3);
    Matrix z = Matrix(3, 1, z_data, 3);
    Matrix g = Matrix(3, 1, g_data, 3);
    Matrix s = Matrix(3, 1, s_data, 3);
    Matrix G = Matrix(3, 3, G_data, 9);
    Matrix P = Matrix(3, 3, P_data, 9);

    Matrix K(3, 3);
    Matrix Y(3, 1);
    Matrix Pout(3, 3);

    MeasUpdate(x, z, g, s, G, P, 3, K, Y, Pout);

    _assert(fabs(K(1, 1) - 0.990099009900990) < TOL_ &&
            fabs(K(2, 2) - 0.990099009900990) < TOL_ &&
            fabs(K(3, 3) - 0.990099009900990) < TOL_);
    _assert(fabs(Y(1) - 1.099009900990099) < TOL_ &&
            fabs(Y(2) - 2.099009900990099) < TOL_ &&
            fabs(Y(3) - 2.900990099009901) < TOL_);
    _assert(fabs(Pout(1, 1) - 0.009900990099010) < TOL_ &&
            fabs(Pout(2, 2) - 0.009900990099010) < TOL_ &&
            fabs(Pout(3, 3) - 0.009900990099010) < TOL_);
    return 0;
}

int TimeUpdate_01() {
    double P_data[36] = {
            100000000, 0, 0, 0, 0, 0,
            0, 100000000, 0, 0, 0, 0,
            0, 0, 100000000, 0, 0, 0,
            0, 0, 0, 1000, 0, 0,
            0, 0, 0, 0, 1000, 0,
            0, 0, 0, 0, 0, 1000
    };
    Matrix P(6, 6, P_data, 36);

    double dt = 1.0;
    double Phi_data[36] = {
            1, 0, 0, dt, 0, 0,
            0, 1, 0, 0, dt, 0,
            0, 0, 1, 0, 0, dt,
            0, 0, 0, 1, 0, 0,
            0, 0, 0, 0, 1, 0,
            0, 0, 0, 0, 0, 1
    };
    Matrix Phi(6, 6, Phi_data, 36);

    double q = 0.01;
    Matrix Qdt(6, 6);
    for (int i = 1; i <= 3; ++i) {
        Qdt(i, i) = 0.25 * q;
        Qdt(i + 3, i + 3) = q;
        Qdt(i, i + 3) = 0.5 * q;
        Qdt(i + 3, i) = 0.5 * q;
    }

    Matrix Pout = TimeUpdate(P, Phi, Qdt);
    _assert(fabs(Pout(1, 1) - 100001000.0025) < TOL_ &&
            fabs(Pout(4, 4) - 1000.01) < TOL_);
    return 0;
}

int sign_01() {
    _assert(fabs(sign_(1e-308, -1e-308) + 1e-308) < TOL_ &&
            fabs(sign_(3, -0.0) - 3) < TOL_);
    return 0;
}

int EqnEquinox_01() {
    _assert(fabs(EqnEquinox(0.0) - 0.0000255262176008994) < TOL_ &&
            fabs(EqnEquinox(58000.0) + 0.0000431583150818001) < TOL_ &&
            fabs(EqnEquinox(2451545.0) + 0.0000670392825380816) < TOL_);
    return 0;
}

int AccelHarmonic_01() {
    double ra[3] = {2000.0, 5500.0, 0.0};
    Matrix r = Matrix(3, 1, ra, 3);
    Matrix E = Matrix::eye(3);

    int n_max = 2;
    int m_max = 2;

    Matrix m = AccelHarmonic(r, E, n_max, m_max);

    _assert(fabs(m(1) * 1e-10 + 0.757459618355740) < TOL_ &&
            fabs(m(2) * 1e-10 + 2.090747514802960) < TOL_ &&
            fabs(m(3) * 1e-10 - 0.000006960610096) < TOL_);
    return 0;
}

int gmst_01() {
    _assert(fabs(gmst(0.0) - 0.973208148169487) < TOL_ &&
            fabs(gmst(58000.0) - 5.9918410478393) < TOL_ &&
            fabs(gmst(2451545.0) - 1.680754521909371) < TOL_ * 10e2);
    return 0;
}

int gast_01() {
    _assert(fabs(gast(0.0) - 0.973233674387088) < TOL_ &&
            fabs(gast(58000.0) - 5.99179788952422) < TOL_ &&
            fabs(gast(-58000.0) - 2.23778075040219) < TOL_ &&
            fabs(gast(2451545.0) - 1.68068748262683) < TOL_ * 100);
    return 0;
}

int G_AccelHarmonic_01() {
    double ra[3] = {7000000.0, 0.0, 0.0};
    Matrix r = Matrix(3, 1, ra, 3);
    Matrix U = Matrix::eye(3);

    int n_max = 2;
    int m_max = 2;

    Matrix G = G_AccelHarmonic(r, U, n_max, m_max);

    _assert(fabs(G(1, 1) - 2.33052264420053e-06) < TOL_ &&
            fabs(G(1, 2) - 2.09290362818138e-11) < TOL_ &&
            fabs(G(1, 3) - 5.32907051820075e-15) < TOL_);
    _assert(fabs(G(2, 1) - 2.0929593988635e-11) < TOL_ &&
            fabs(G(2, 2) + 1.16369909716372e-06) < TOL_ &&
            fabs(G(2, 3) - 5.47310677681545e-15) < TOL_);
    _assert(fabs(G(3, 1) - 3.34004807862169e-15) < TOL_ &&
            fabs(G(3, 2) - 5.47310880258077e-15) < TOL_ &&
            fabs(G(3, 3) + 1.16682354671017e-06) < TOL_);
    return 0;
}

int LTC_01() {
    // Kaena Point station
    double lon = Rad * (-158.2706);
    double lat = Rad * 21.5748;

    Matrix M = LTC(lon, lat);

    _assert(fabs(M(1, 1) - 0.370223471399199) < TOL_ &&
            fabs(M(1, 2) + 0.928942722252092) < TOL_ &&
            fabs(M(1, 3) - 0.0) < TOL_);
    _assert(fabs(M(2, 1) - 0.341586711932422) < TOL_ &&
            fabs(M(2, 2) - 0.136136938528208) < TOL_ &&
            fabs(M(2, 3) - 0.929938305587722) < TOL_);
    _assert(fabs(M(3, 1) + 0.863859421119156) < TOL_ &&
            fabs(M(3, 2) + 0.344284987681776) < TOL_ &&
            fabs(M(3, 3) - 0.367715580035218) < TOL_);
    return 0;
}

int NutMatrix_01() {
    Matrix NutMat = NutMatrix(MJD_J2000 + 1.0 / 86400);

    _assert(fabs(NutMat(1, 1) - 0.999999997721708) < TOL_ &&
            fabs(NutMat(1, 2) - 6.19323105974282e-05) < TOL_ &&
            fabs(NutMat(1, 3) - 2.68509428011841e-05) < TOL_);
    _assert(fabs(NutMat(2, 1) + 6.19330621904886e-05) < TOL_ &&
            fabs(NutMat(2, 2) - 0.999999997690389) < TOL_ &&
            fabs(NutMat(2, 3) - 2.79913820497302e-05) < TOL_);
    _assert(fabs(NutMat(3, 1) + 2.68492091682017e-05) < TOL_ &&
            fabs(NutMat(3, 2) + 2.79930449471055e-05) < TOL_ &&
            fabs(NutMat(3, 3) - 0.999999999247755) < TOL_);
    return 0;
}

int PoleMatrix_01() {
    Matrix PoleMat = PoleMatrix(-5.59378710288119e-07, 2.33559844213816e-06);

    _assert(fabs(PoleMat(1, 1) - 0.999999999999844) < TOL_ &&
            fabs(PoleMat(1, 2) + 1.30648404431293e-12) < TOL_ &&
            fabs(PoleMat(1, 3) + 5.59378710286564e-07) < TOL_);
    _assert(fabs(PoleMat(2, 1) - 0.0) < TOL_ &&
            fabs(PoleMat(2, 2) - 0.999999999997273) < TOL_ &&
            fabs(PoleMat(2, 3) + 2.33559844213604e-06) < TOL_);
    _assert(fabs(PoleMat(3, 1) - 5.5937871028809e-07) < TOL_ &&
            fabs(PoleMat(3, 2) - 2.33559844213567e-06) < TOL_ &&
            fabs(PoleMat(3, 3) - 0.999999999997116) < TOL_);
    return 0;
}

int PrecMatrix_01() {
    Matrix PrecMat = PrecMatrix(MJD_J2000, 49746.1163657406);

    _assert(fabs(PrecMat(1, 1) - 0.999999279432408) < TOL_ &&
            fabs(PrecMat(1, 2) - 0.00110100876673246) < TOL_ &&
            fabs(PrecMat(1, 3) - 0.000478449956744648) < TOL_);
    _assert(fabs(PrecMat(2, 1) + 0.00110100876673353) < TOL_ &&
            fabs(PrecMat(2, 2) - 0.999999393889629) < TOL_ &&
            fabs(PrecMat(2, 3) + 2.63386664233359e-07) < TOL_);
    _assert(fabs(PrecMat(3, 1) + 0.000478449956742193) < TOL_ &&
            fabs(PrecMat(3, 2) + 2.6339112237367e-07) < TOL_ &&
            fabs(PrecMat(3, 3) - 0.999999885542778) < TOL_);
    return 0;
}

int GHAMatrix_01() {
    Matrix GHAMat = GHAMatrix(49746.1163657406);

    _assert(fabs(GHAMat(1, 1) + 0.9841444735756) < TOL_ &&
            fabs(GHAMat(1, 2) - 0.177368698282998) < TOL_ &&
            fabs(GHAMat(1, 3) - 0.0) < TOL_);
    _assert(fabs(GHAMat(2, 1) + 0.177368698282998) < TOL_ &&
            fabs(GHAMat(2, 2) + 0.9841444735756) < TOL_ &&
            fabs(GHAMat(2, 3) - 0.0) < TOL_);
    _assert(fabs(GHAMat(3, 1) - 0.0) < TOL_ &&
            fabs(GHAMat(3, 2) - 0.0) < TOL_ &&
            fabs(GHAMat(3, 3) - 1.0) < TOL_);
    return 0;
}

int Geodetic_01() {
    double lon, lat, h;
    double a[3] = {1.0, 2.0, 3.0};
    Matrix r = Matrix(3, 1, a, 3);
    Geodetic(r, lon, lat, h);

    _assert(fabs(lon - 1.10714871779409) < TOL_ &&
            fabs(lat - 1.57074413624392) < TOL_ &&
            fabs(h * 1e-6 + 6.3567486165338) < TOL_);
    return 0;
}

int IERS_01() {

    double x_pole, y_pole, UT1_UTC, LOD, dpsi, deps, dx_pole, dy_pole, TAI_UTC;

    IERS(eop, 37670, x_pole, y_pole, UT1_UTC, LOD, dpsi, deps, dx_pole, dy_pole, TAI_UTC, 'l');

    _assert(fabs(x_pole + 1.33803727849421e-07) < TOL_ &&
            fabs(y_pole - 1.05835311399893e-06) < TOL_ &&
            fabs(UT1_UTC - 0.030535300000000) < TOL_ &&
            fabs(LOD - 0.001382000000000) < TOL_ &&
            fabs(dpsi - 3.094371801049724e-07) < TOL_ &&
            fabs(deps - 2.930698702307145e-08) < TOL_ &&
            fabs(dx_pole - 0.0) < TOL_ &&
            fabs(dy_pole - 0.0) < TOL_ &&
            fabs(TAI_UTC - 2.0) < TOL_);
    return 0;
}

int JPL_Eph_DE430_01() {
    Matrix r_Mercury(3, 1), r_Venus(3, 1), r_Earth(3, 1), r_Mars(3, 1),
            r_Jupiter(3, 1), r_Saturn(3, 1), r_Uranus(3, 1), r_Neptune(3, 1),
            r_Pluto(3, 1), r_Moon(3, 1), r_Sun(3, 1);
    JPL_Eph_DE430(49746.1107720813, r_Mercury, r_Venus, r_Earth, r_Mars,
                  r_Jupiter, r_Saturn, r_Uranus, r_Neptune, r_Pluto,
                  r_Moon, r_Sun);

    _assert(fabs(r_Mercury(1) * 1e-10 - 8.377907564855100) < TOL_ &&
            fabs(r_Mercury(2) * 1e-10 + 6.529205179483268) < TOL_ &&
            fabs(r_Mercury(3) * 1e-10 + 2.339225533913997) < TOL_);
    _assert(fabs(r_Venus(1) * 1e-11 + 0.152322317037660) < TOL_ &&
            fabs(r_Venus(2) * 1e-11 + 1.101334280274653) < TOL_ &&
            fabs(r_Venus(3) * 1e-11 + 0.410210661150587) < TOL_);
    _assert(fabs(r_Earth(1) * 1e-11 + 0.924684582375320) < TOL_ &&
            fabs(r_Earth(2) * 1e-11 - 1.063967349682521) < TOL_ &&
            fabs(r_Earth(3) * 1e-11 - 0.461309275429331) < TOL_);
    _assert(fabs(r_Mars(1) * 1e-10 + 8.827926215704085) < TOL_ &&
            fabs(r_Mars(2) * 1e-10 - 4.696446520285535) < TOL_ &&
            fabs(r_Mars(3) * 1e-10 - 2.907088764515343) < TOL_);
    _assert(fabs(r_Jupiter(1) * 1e-11 + 2.983896254699393) < TOL_ &&
            fabs(r_Jupiter(2) * 1e-11 + 7.544995288385803) < TOL_ &&
            fabs(r_Jupiter(3) * 1e-11 + 3.144110429896569) < TOL_);
    _assert(fabs(r_Saturn(1) * 1e-12 - 1.482031269946539) < TOL_ &&
            fabs(r_Saturn(2) * 1e-12 + 0.453875617598394) < TOL_ &&
            fabs(r_Saturn(3) * 1e-12 + 0.249403400202920) < TOL_);
    _assert(fabs(r_Uranus(1) * 1e-12 - 1.412364844217200) < TOL_ &&
            fabs(r_Uranus(2) * 1e-12 + 2.511357129966920) < TOL_ &&
            fabs(r_Uranus(3) * 1e-12 + 1.118109547729565) < TOL_);
    _assert(fabs(r_Neptune(1) * 1e-12 - 1.871247743893495) < TOL_ &&
            fabs(r_Neptune(2) * 1e-12 + 3.928978347157190) < TOL_ &&
            fabs(r_Neptune(3) * 1e-12 + 1.655021340136643) < TOL_);
    _assert(fabs(r_Pluto(1) * 1e-12 + 2.171417806084585) < TOL_ &&
            fabs(r_Pluto(2) * 1e-12 + 3.915434641172163) < TOL_ &&
            fabs(r_Pluto(3) * 1e-12 + 0.552716789894166) < TOL_);
    _assert(fabs(r_Moon(1) * 1e-08 - 0.892732558617431) < TOL_ &&
            fabs(r_Moon(2) * 1e-08 + 3.366262378225644) < TOL_ &&
            fabs(r_Moon(3) * 1e-08 + 1.146634306237374) < TOL_);
    _assert(fabs(r_Sun(1) * 1e-11 - 0.922957499797460) < TOL_ &&
            fabs(r_Sun(2) * 1e-11 + 1.053770129120909) < TOL_ &&
            fabs(r_Sun(3) * 1e-11 + 0.456871550046014) < TOL_);
    return 0;
}

int Accel_01() {
    auxParam.Mjd_UTC = 49746.1163541665;
    auxParam.n = 20;
    auxParam.m = 20;
    auxParam.sun = true;
    auxParam.moon = true;
    auxParam.planets = true;
    double Yarray[] = {6221397.62857869, 2867713.77965738,
                       3006155.98509949, 4645.04725161807,
                       -2752.21591588205, -7507.99940987033};
    Matrix Y = Matrix(6, 1, Yarray, 6);

    Matrix AccelMat = Accel(MJD_J2000, Y);

    _assert(fabs(AccelMat(1) - 4645.04725161807) < TOL_ &&
            fabs(AccelMat(2) + 2752.21591588205) < TOL_ &&
            fabs(AccelMat(3) + 7507.99940987033) < TOL_);
    _assert(fabs(AccelMat(4) + 5.92429522874995) < TOL_ &&
            fabs(AccelMat(5) + 2.73078680171038) < TOL_ &&
            fabs(AccelMat(6) + 2.86933058377641) < TOL_);
    return 0;
}

int Accel_02() {
    auxParam.Mjd_UTC = 49746.1163541665;
    auxParam.Mjd_TT = 49746.1170623147;
    auxParam.n = 20;
    auxParam.m = 20;
    auxParam.sun = true;
    auxParam.moon = true;
    auxParam.planets = true;

    double Yarray[] = {6878000.0, 1000.0, -2000.0, -200.0, 7550.0, 100.0};
    Matrix Y(6, 1, Yarray, 6);
    Matrix AccelMat = Accel(MJD_J2000, Y);

    _assert(fabs(AccelMat(1) + 200.0) < TOL_ &&
            fabs(AccelMat(2) - 7550) < TOL_ &&
            fabs(AccelMat(3) - 100.0) < TOL_);
    _assert(fabs(AccelMat(4) + 8.43765658331181) < TOL_ &&
            fabs(AccelMat(5) + 0.00115924469366841) < TOL_ &&
            fabs(AccelMat(6) - 0.00251874989793148) < TOL_);
    return 0;
}

int Accel_03() {
    auxParam.Mjd_UTC = 53000.0;
    auxParam.n = 10;
    auxParam.m = 10;
    auxParam.sun = false;
    auxParam.moon = true;
    auxParam.planets = false;

    double Yarray[] = {6878000.0, 1000.0, -2000.0, -200.0, 7550.0, 100.0};
    Matrix Y(6, 1, Yarray, 6);
    Matrix AccelMat = Accel(auxParam.Mjd_UTC, Y);

    _assert(fabs(AccelMat(1) + 200.0) < TOL_ &&
            fabs(AccelMat(2) - 7550) < TOL_ &&
            fabs(AccelMat(3) - 100.0) < TOL_ &&
            fabs(AccelMat(4) + 8.43750708231707) < TOL_ &&
            fabs(AccelMat(5) + 0.00131889962222729) < TOL_ &&
            fabs(AccelMat(6) - 0.00247717103915464) < TOL_);

    return 0;
}


int VarEqn_01() {
    auxParam.Mjd_UTC = 58000;
    auxParam.Mjd_TT = 58000.1;
    auxParam.n = 20;
    auxParam.m = 20;
    auxParam.sun = true;
    auxParam.moon = true;
    auxParam.planets = true;
    double Yarray[] = {7081724.51601926, 1384037.1930286,
                       208050.100775178, 794.267807462975,
                       -3043.46324226258, -6732.60129776107,
                       1.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                       0.0, 1.0, 0.0, 0.0, 0.0, 0.0,
                       0.0, 0.0, 1.0, 0.0, 0.0, 0.0,
                       0.0, 0.0, 0.0, 1.0, 0.0, 0.0,
                       0.0, 0.0, 0.0, 0.0, 1.0, 0.0,
                       0.0, 0.0, 0.0, 0.0, 0.0, 1.0
    };
    Matrix Y = Matrix(42, 1, Yarray, 42);
    Matrix VarEqnMat = VarEqn(MJD_J2000, Y);

    _assert(fabs(VarEqnMat(1) - 794.267807462975) < TOL_ &&
            fabs(VarEqnMat(2) + 3043.46324226258) < TOL_ &&
            fabs(VarEqnMat(3) + 6732.60129776107) < TOL_);
    _assert(fabs(VarEqnMat(4) + 7.51367051488493) < TOL_ &&
            fabs(VarEqnMat(5) + 1.4684498246214) < TOL_ &&
            fabs(VarEqnMat(6) + 0.221298652184411) < TOL_);
    _assert(fabs(VarEqnMat(7) - 0.0) < TOL_ &&
            fabs(VarEqnMat(8) - 0.0) < TOL_ &&
            fabs(VarEqnMat(9) - 0.0) < TOL_);
    _assert(fabs(VarEqnMat(10) - 2.00494108870686e-06) < TOL_ &&
            fabs(VarEqnMat(11) - 5.99213035634705e-07) < TOL_ &&
            fabs(VarEqnMat(12) - 9.04438409721209e-08) < TOL_);
    _assert(fabs(VarEqnMat(13) - 0.0) < TOL_ &&
            fabs(VarEqnMat(14) - 0.0) < TOL_ &&
            fabs(VarEqnMat(15) - 0.0) < TOL_);
    _assert(fabs(VarEqnMat(16) - 5.99213043628311e-07) < TOL_ &&
            fabs(VarEqnMat(17) + 9.43921478091525e-07) < TOL_ &&
            fabs(VarEqnMat(18) - 1.76749646862984e-08) < TOL_);
    _assert(fabs(VarEqnMat(19) - 0.0) < TOL_ &&
            fabs(VarEqnMat(20) - 0.0) < TOL_ &&
            fabs(VarEqnMat(21) - 0.0) < TOL_);
    _assert(fabs(VarEqnMat(22) - 9.04438435256338e-08) < TOL_ &&
            fabs(VarEqnMat(23) - 1.76749652691655e-08) < TOL_ &&
            fabs(VarEqnMat(24) + 1.06101961272476e-06) < TOL_);
    _assert(fabs(VarEqnMat(25) - 1.0) < TOL_ &&
            fabs(VarEqnMat(26) - 0.0) < TOL_ &&
            fabs(VarEqnMat(27) - 0.0) < TOL_);
    _assert(fabs(VarEqnMat(28) - 0.0) < TOL_ &&
            fabs(VarEqnMat(29) - 0.0) < TOL_ &&
            fabs(VarEqnMat(30) - 0.0) < TOL_);
    _assert(fabs(VarEqnMat(31) - 0.0) < TOL_ &&
            fabs(VarEqnMat(32) - 1.0) < TOL_ &&
            fabs(VarEqnMat(33) - 0.0) < TOL_);
    _assert(fabs(VarEqnMat(34) - 0.0) < TOL_ &&
            fabs(VarEqnMat(35) - 0.0) < TOL_ &&
            fabs(VarEqnMat(36) - 0.0) < TOL_);
    _assert(fabs(VarEqnMat(37) - 0.0) < TOL_ &&
            fabs(VarEqnMat(38) - 0.0) < TOL_ &&
            fabs(VarEqnMat(39) - 1.0) < TOL_);
    _assert(fabs(VarEqnMat(40) - 0.0) < TOL_ &&
            fabs(VarEqnMat(41) - 0.0) < TOL_ &&
            fabs(VarEqnMat(42) - 0.0) < TOL_);
    return 0;
}

int DEInteg_01() {
    auxParam.Mjd_UTC = 49746.1112847221;
    auxParam.n = 20;
    auxParam.m = 20;
    auxParam.sun = true;
    auxParam.moon = true;
    auxParam.planets = true;

    double Y0_apra[] = {6221397.62857869,
                        2867713.77965738,
                        3006155.98509949,
                        4645.04725161807,
                        -2752.21591588205,
                        -7507.99940987033};
    Matrix Y0_apr = Matrix(6, 1, Y0_apra, 6);

    DEInteg(Accel, 0, -134.999991953373, 1e-13, 1e-6, 6, Y0_apr);

    _assert(fabs(Y0_apr(1) - 5542555.93722861) < TOL_ * 10e6);
    _assert(fabs(Y0_apr(2) - 3213514.8673492) < TOL_ * 10e6);
    _assert(fabs(Y0_apr(3) - 3990892.97587686) < TOL_ * 10e6);
    _assert(fabs(Y0_apr(4) - 5394.06842166353) < TOL_ * 10e3);
    _assert(fabs(Y0_apr(5) + 2365.21337882342) < TOL_ * 10e3);
    _assert(fabs(Y0_apr(6) + 7061.84554200298) < TOL_ * 10e3);

    return 0;
}

int DEInteg_02() {
    auxParam.Mjd_UTC = 49746.1163541665;

    double Y_apra[] = {
            7101576.98990384,
            1295199.87127754,
            12739.2823333892,
            576.004651193009,
            -3084.62203617271,
            -6736.02594582756,
            1.00002525535511,
            7.08259815373561e-06,
            1.91608861002907e-07,
            1.01043851887223e-05,
            2.82768336557965e-06,
            6.44131451075285e-08,
            7.08259834024473e-06,
            0.999988040046622,
            3.53015288644891e-08,
            2.82768357826951e-06,
            -4.78603729288896e-06,
            1.18527461137171e-08,
            1.9160935046062e-07,
            3.53016114843062e-08,
            0.999986704774626,
            6.44136325079115e-08,
            1.18528331537947e-08,
            -5.31820682446032e-06,
            5.00001498082565,
            1.1781862826826e-05,
            2.68389762645616e-07,
            1.00002526606744,
            7.05571100144785e-06,
            1.30455137405173e-07,
            1.17818628919961e-05,
            4.99995293819715,
            4.93630678596596e-08,
            7.05571117883108e-06,
            0.999988029832331,
            2.39618837211068e-08,
            2.68390168073246e-07,
            4.93631303180711e-08,
            4.99995072081276,
            1.30455621823661e-07,
            2.39619698989173e-08,
            0.999986704276552
    };
    Matrix Y0_apr = Matrix(42, 1, Y_apra, 42);

    DEInteg(VarEqn, 0, 4.99997287988663, 1e-13, 1e-6, 42, Y0_apr);

    return 0;
}

int angl_01(){
    double vec1array[] = {1.0,-5.5,2.0};
    double vec2array[] = {-2.2,0.0,7.5};
    Matrix vec1(3,1, vec1array, 3);
    Matrix vec2(3,1, vec2array, 3);

    _assert(fabs(angl(vec1, vec2) - 1.29134069114822) < TOL_);
    return 0;
}

int gibbs_01(){
    double vec1array[] = {1.0,-5.5,2.0};
    double vec2array[] = {-2.2,0.0,7.5};
    double vec3array[] = {-3.2,5.0,-9.5};
    Matrix vec1(3,1, vec1array, 3);
    Matrix vec2(3,1, vec2array, 3);
    Matrix vec3(3,1, vec3array, 3);
    double theta, theta1, copa;
    Matrix v2(3,1);
    char * error;

    gibbs(vec1, vec2, vec3, v2, theta, theta1, copa, error);

    _assert(fabs(v2(1)*1e-06 + 3.42915690467666) < TOL_ &&
            fabs(v2(2)*1e-06 - 6.62205674314325) < TOL_ &&
            fabs(v2(3)*1e-06 - 1.73969395729477) < TOL_);
    _assert(fabs(theta - 1.29134069114822) < TOL_);
    _assert(fabs(theta1 - 2.39403716882476) < TOL_);
    _assert(fabs(copa - 0.559073570331478) < TOL_);

    return 0;
}

int hgibbs_01(){
    double vec1array[] = {1.0,-5.5,2.0};
    double vec2array[] = {-2.2,0.0,7.5};
    double vec3array[] = {-3.2,5.0,-9.5};
    Matrix vec1(3,1, vec1array, 3);
    Matrix vec2(3,1, vec2array, 3);
    Matrix vec3(3,1, vec3array, 3);
    double theta, theta1, copa;
    Matrix v2(3,1);
    char * error;

    hgibbs(vec1, vec2, vec3, 51000.5,53000.1,49500.5, v2, theta, theta1, copa, error);

    _assert(fabs(v2(1)*1e-20 - 1.07643933612008) < TOL_ &&
            fabs(v2(2)*1e-20 + 2.43533283328675) < TOL_ &&
            fabs(v2(3)*1e-20 + 1.90701170075358) < TOL_);
    _assert(fabs(theta - 1.29134069114822) < TOL_);
    _assert(fabs(theta1 - 2.39403716882476) < TOL_);
    _assert(fabs(copa - 0.559073570331478) < TOL_);

    return 0;
}

int elements_01(){
    double y[6] = {1.0,3.0,-5.5,0.0,3.5,6.4};
    double p, a, e, i, Omega, omega, M;

    elements(y, p, a, e, i, Omega, omega, M);

    _assert(fabs(p - 3.84247573218198e-12) < TOL_);
    _assert(fabs(a - 3.17214438511372) < TOL_);
    _assert(fabs(e - 0.999999999999394) < TOL_);
    _assert(fabs(i - 1.48124454716499) < TOL_);
    _assert(fabs(Omega - 1.40585853086175) < TOL_);
    _assert(fabs(omega - 2.08555969410608) < TOL_);
    _assert(fabs(M - 3.14159404284359) < TOL_);

    return 0;
}

int anglesg_01(){
    double lat = Rad * 21.5748;     // [rad]
    double lon = Rad * (-158.2706); // [rad]
    double alt = 300.20;            // [m]
    Matrix Rs = Position(lon, lat, alt);
    double Mjd1 = obs(1, 1);
    double Mjd2 = obs(9, 1);
    double Mjd3 = obs(18, 1);
    Matrix r2(3, 1);
    Matrix v2(3, 1);

    anglesg(obs(1,2),obs(9,2),obs(18,2),obs(1,3),obs(9,3),obs(18,3), Mjd1,Mjd2,Mjd3,Rs,Rs,Rs,r2,v2);

    r2.print();
    v2.print();

    _assert(fabs(r2(1)*1e-06 - 6.22139762857869) < TOL_ &&
            fabs(r2(2)*1e-06 - 2.86771377965738) < TOL_ &&
            fabs(r2(3)*1e-06 - 3.00615598509949) < TOL_);
    _assert(fabs(v2(1)*1e-03 - 4.64504725161807) < TOL_*10 &&
            fabs(v2(2)*1e-03 + 2.75221591588205) < TOL_*10 &&
            fabs(v2(3)*1e-03 + 7.50799940987033) < TOL_*10);

    return 0;
}

int all_tests() {
    try{
        _verify(Mjday_01);
        _verify(Mjday_02);
        _verify(R_x_01);
        _verify(R_y_01);
        _verify(R_z_01);
        _verify(Legendre_01);
        _verify(Frac_01);
        _verify(Position_01);
        _verify(AccelPointMass_01);
        _verify(AccelPointMass_02);
        _verify(Mjd_TDB_01);
        _verify(Mjd_TDB_02);
        _verify(NutAngles_01);
        _verify(NutAngles_02);
        _verify(MeanObliquity_01);
        _verify(MeanObliquity_02);
        _verify(Cheb3D_01);
        _verify(Cheb3D_02);
        _verify(timediff_01);
        _verify(timediff_02);
        _verify(AzElPa_01);
        _verify(unit_01);
        _verify(unit_02);
        _verify(MeasUpdate_01);
        _verify(TimeUpdate_01);
        _verify(sign_01);
        _verify(EqnEquinox_01);
        _verify(AccelHarmonic_01);
        _verify(gmst_01);
        _verify(gast_01);
        _verify(G_AccelHarmonic_01);
        _verify(LTC_01);
        _verify(NutMatrix_01);
        _verify(PoleMatrix_01);
        _verify(PrecMatrix_01);
        _verify(Geodetic_01);
        _verify(IERS_01);
        _verify(GHAMatrix_01);
        _verify(JPL_Eph_DE430_01);
        _verify(Accel_01);
        _verify(Accel_02);
        _verify(Accel_03);
        _verify(VarEqn_01);
        _verify(DEInteg_01);
        _verify(DEInteg_02);
        _verify(angl_01);
        _verify(gibbs_01);
        _verify(hgibbs_01);
        _verify(elements_01);
        _verify(anglesg_01);
    } catch (exception& e) {
        cout << e.what() << endl;
        FAIL();
        return 1;
    }

    return 0;
}

int main() {
    auxParam.Mjd_UTC = 49746.1163541665;
    auxParam.n = 20;
    auxParam.m = 20;
    auxParam.sun = true;
    auxParam.moon = true;
    auxParam.planets = true;
    auxParam.Mjd_TT = 49746.1170623147;

    loadEOP("./data/eop19620101.txt");
    loadDE430Coeff("./data/DE430Coeff.txt");
    loadGMS03S("./data/GGM03S.txt");
    loadGEOS3("data/GEOS3.txt");

    int result = all_tests();

    if (result == 0)
        printf("\nPASSED\n");
    else
        FAIL();

    printf("Tests run: %d\n", tests_run);

    return result != 0;
}