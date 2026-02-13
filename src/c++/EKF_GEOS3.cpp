//--------------------------------------------------------------------------
//
// Initial Orbit Determination using Gauss and Extended Kalman Filter methods
//
// References:
//   O. Montenbruck, E. Gill, "Satellite Orbits - Models, Methods, and
//   Applications", Springer Verlag, Heidelberg, 2000.
//
//   D. Vallado, "Fundamentals of Astrodynamics and Applications",
//   4th Edition, 2013.
//
//   G. Seeber, "Satellite Geodesy", 2nd Edition, 2003.
//
// Last modified:   2020/03/16   Meysam Mahooti
//--------------------------------------------------------------------------

// Note: This is a migration of the original MatLab code to C++, by Martín Hernández Tonzán

#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iomanip>
#include <chrono>
#include "include/SAT_Const.h"
#include "include/Mjday.h"
#include "include/Matrix.h"
#include "include/Position.h"
#include "include/Accel.h"
#include "include/VarEqn.h"
#include "include/DEInteg.h"
#include "include/LTC.h"
#include "include/IERS.h"
#include "include/timediff.h"
#include "include/gmst.h"
#include "include/TimeUpdate.h"
#include "include/R_z.h"
#include "include/AzElPa.h"
#include "include/MeasUpdate.h"
#include "EKF_Global.h"
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

int main() {

    int nobs = 46;

    auto start = std::chrono::high_resolution_clock::now();

    loadEOP("./data/eop19620101.txt");
    loadDE430Coeff("./data/DE430Coeff.txt");
    loadGMS03S("./data/GGM03S.txt");
    loadGEOS3("./data/GEOS3.txt");

    auto endficheros = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> durationficheros = endficheros - start;
    std::cout << "Tiempo de carga de ficheros: " << durationficheros.count() / 1000 << " s" << std::endl;


    double sigma_range = 92.5;      // [m]
    double sigma_az = 0.0224 * Rad; // [rad]
    double sigma_el = 0.0139 * Rad; // [rad]

    // Kaena Point station
    double lat = Rad * 21.5748;     // [rad]
    double lon = Rad * (-158.2706); // [rad]
    double alt = 300.20;                // [m]

    Matrix Rs = Position(lon, lat, alt);

    double Mjd1 = obs(1, 1);
    //double Mjd1 = obs[0][0];
    double Mjd2 = obs(9, 1);
    //double Mjd2 = obs[8][0];
    double Mjd3 = obs(18, 1);
    //double Mjd3 = obs[17][0];

    Matrix r2(3,1); Matrix v2(3,1);

    anglesg(obs(1,2),obs(9,2),obs(18,2),obs(1,3),obs(9,3),obs(18,3),Mjd1,Mjd2,Mjd3,Rs,Rs,Rs, r2, v2);

    //anglesg(obs[0][1],obs[8][1],obs[17][1],obs[0][2],obs[8][2],obs[17][2],Mjd1,Mjd2,Mjd3,Rs,Rs,Rs, r2, v2);
    //anglesdr(obs(1,2),obs(9,2),obs(18,2),obs(1,3),obs(9,3),obs(18,3),...
    //                    Mjd1,Mjd2,Mjd3,Rs,Rs,Rs,r2,v2);

    //Y0_apr = [r2;v2];

    // profe
    /*double Y0_aprArray[] = {6173608.86411755, 2629396.63515668, 2920149.6667258,
                            4080.5957475354, -4835.33102038702, -8123.49975545833};*/

    double Y0_aprArray[6] = {r2(1), r2(2), r2(3), v2(1), v2(2), v2(3)};

    // normal
   // double Y0_aprArray[] = {6221397.62857869, 2867713.77965738, 3006155.98509949,
     //                       4645.04725161801, -2752.21591588201, -7507.99940987021};

    double Mjd0 = Mjday(1995, 1, 29, 02, 38, 0);

    double Mjd_UTC = obs(9, 1);
    //double Mjd_UTC = obs[8][0];

    auxParam.Mjd_UTC = Mjd_UTC;
    auxParam.n = 20;
    auxParam.m = 20;
    auxParam.sun = true;
    auxParam.moon = true;
    auxParam.planets = true;

    int n_eqn = 6;

    Matrix Y0_apr = Matrix(6, 1, Y0_aprArray, 6);
    Y0_apr.print();

    DEInteg(Accel, 0, -(obs(9, 1) - Mjd0) * 86400.0, 1e-13, 1e-6, 6, Y0_apr);
    //DEInteg(Accel, 0, -(obs[8][0] - Mjd0) * 86400.0, 1e-13, 1e-6, 6, Y0_apr);
    Matrix Y = Y0_apr;

    Matrix P = Matrix(6, 6);

    for (int i = 1; i <= 3; i++) {
        P(i, i) = 1e8;
    }
    for (int i = 4; i <= 6; i++) {
        P(i, i) = 1e3;
    }

    Matrix LT = LTC(lon, lat);

    Matrix yPhi(42, 1);
    Matrix Phi(6, 6);

    // Measurement loop
    double t = 0.0;
    double x_pole, y_pole, UT1_UTC, LOD, dpsi, deps, dx_pole, dy_pole, TAI_UTC;
    double UT1_TAI, UTC_GPS, UT1_GPS, TT_UTC, GPS_UTC;
    double Azim, Elev;
    Matrix dAds(3, 1), dEds(3, 1), K(3, 3);

    for (int i = 1; i <= nobs; i++) {
        // Previous step
        double t_old = t;
        Matrix Y_old = Y;

        // Time increment and propagation
        Mjd_UTC = obs(i, 1);                  // Modified Julian Date
        //Mjd_UTC = obs[i][0];                  // Modified Julian Date
        t = (Mjd_UTC - Mjd0) * 86400.0;         // Time since epoch [s]

        IERS(eop,Mjd_UTC, x_pole, y_pole, UT1_UTC, LOD, dpsi, deps, dx_pole, dy_pole, TAI_UTC, 'l');
        timediff(UT1_UTC, TAI_UTC, UT1_TAI, UTC_GPS, UT1_GPS, TT_UTC, GPS_UTC);
        double Mjd_TT = Mjd_UTC + TT_UTC / 86400;
        double Mjd_UT1 = Mjd_TT + (UT1_UTC - TT_UTC) / 86400.0;
        auxParam.Mjd_UTC = Mjd_UTC;
        auxParam.Mjd_TT = Mjd_TT;

        for (int ii = 1; ii <= 6; ii++) {
            yPhi(ii) = Y_old(ii);
            for (int j = 1; j <= 6; j++) {
                if (ii == j)
                    yPhi(6 * j + ii) = 1.0;
                else
                    yPhi(6 * j + ii) = 0.0;
            }
        }

        DEInteg(VarEqn, 0, t - t_old, 1e-13, 1e-6, 42, yPhi);

        // Extract state transition matrices
        for (int j = 1; j <= 6; j++) {
            for (int ix = 1; ix <= 6; ++ix) {
                Phi(ix, j) = yPhi(6 * j + ix);
            }
        }

        DEInteg(Accel, 0, t - t_old, 1e-13, 1e-6, 6, Y_old);
        Y = Y_old;

        // Topocentric coordinates
        double theta = gmst(Mjd_UT1);                    // Earth rotation
        Matrix U = R_z(theta);
        Matrix r = Y.subvector(1, 3).transpose();
        Matrix s = LT * (U * r - Rs);                          // Topocentric position [m]

        // Time update
        P = TimeUpdate(P, Phi);

        // Azimuth and partials
        AzElPa(s, Azim, Elev, dAds, dEds);     // Azimuth, Elevation
        Matrix dAdY = (dAds.transpose() * LT * U).concatRow(Matrix(1, 3));

        Matrix obsAux(1, 1);
        obsAux(1) = obs(i, 2);
        //obsAux(1) = obs[i][1];
        Matrix AzimAux(1, 1);
        AzimAux(1) = Azim;
        Matrix SigmaAux(1, 1);
        SigmaAux(1) = sigma_az;

        // Measurement update
        MeasUpdate(Y, obsAux, AzimAux, SigmaAux, dAdY, P, 6, K, Y, P);

        // Elevation and partials
        r = Y.subvector(1, 3).transpose();
        s = LT * (U * r - Rs);                              // Topocentric position [m]
        AzElPa(s, Azim, Elev, dAds, dEds);     // Azimuth, Elevation
        Matrix dEdY = (dEds.transpose() * LT * U).concatRow(Matrix(1, 3));

        obsAux(1) = obs(i, 3);
        //obsAux(1) = obs[i][2];
        Matrix ElevAux(1, 1);
        ElevAux(1) = Elev;
        SigmaAux(1) = sigma_el;

        // Measurement update
        MeasUpdate(Y, obsAux, ElevAux, SigmaAux, dEdY, P, 6, K, Y, P);

        // Range and partials
        r = Y.subvector(1, 3).transpose();
        s = LT * (U * r - Rs);                          // Topocentric position [m]
        double Dist = Matrix::norm(s);
        Matrix dDds = (s / Dist).transpose();         // Range
        Matrix dDdY = (dDds * LT * U).concatRow(Matrix(1, 3));

        obsAux(1) = obs(i, 4);
        //obsAux(1) = obs[i][3];
        Matrix DistAux(1, 1);
        DistAux(1) = Dist;
        SigmaAux(1) = sigma_range;

        // Measurement update
        MeasUpdate(Y, obsAux, DistAux, SigmaAux, dDdY, P, 6, K, Y, P);
    }

    IERS(eop, obs(46, 1), x_pole, y_pole, UT1_UTC, LOD, dpsi, deps, dx_pole, dy_pole, TAI_UTC, 'l');
    //IERS(obs[45][0], x_pole, y_pole, UT1_UTC, LOD, dpsi, deps, dx_pole, dy_pole, TAI_UTC, 'l');
    timediff(UT1_UTC, TAI_UTC, UT1_TAI, UTC_GPS, UT1_GPS, TT_UTC, GPS_UTC);
    double Mjd_TT = Mjd_UTC + TT_UTC / 86400;
    auxParam.Mjd_UTC = Mjd_UTC;
    auxParam.Mjd_TT = Mjd_TT;
    DEInteg(Accel, 0, -(obs(46, 1) - obs(1, 1)) * 86400.0, 1e-13, 1e-6, 6, Y);
    // DEInteg(Accel, 0, -(obs[45][0] - obs[0][0]) * 86400.0, 1e-13, 1e-6, 6, Y);

    Matrix Y0 = Y;

    double Y_trueArray[] = {5753.173e3, 2673.361e3, 3440.304e3, 4.324207e3, -1.924299e3, -5.728216e3};

    Matrix Y_true(6, 1, Y_trueArray, 6);

    cout << setprecision(10) << "Error of Position Estimation";
    cout << endl << "dX\t" << Y0(1) - Y_true(1) << " [m]";
    cout << endl << "dY\t" << Y0(2) - Y_true(2) << " [m]";
    cout << endl << "dz\t" << Y0(3) - Y_true(3) << " [m]";
    cout << setprecision(10) << endl << "Error of Velocity Estimation";
    cout << endl << "dVx\t" << Y0(4) - Y_true(4) << " [m/s]";
    cout << endl << "dVy\t" << Y0(5) - Y_true(5) << " [m/s]";
    cout << endl << "dVz\t" << Y0(6) - Y_true(6) << " [m/s]";

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> duration = end - start;
    std::cout << endl << "Tiempo de ejecución total: " << duration.count() / 1000 << " s" << std::endl;


    return 0;

}