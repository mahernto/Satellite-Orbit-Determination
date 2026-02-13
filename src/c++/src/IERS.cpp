/// @file IERS.cpp
/// @brief Source file for IERS function
/// @date 2025-04-16
/// @author Martín Hernández Tonzán

#include <cmath>
#include <iostream>
#include "../include/IERS.h"

void IERS(const Matrix& eop, double Mjd_UTC, double& x_pole, double& y_pole, double& UT1_UTC, double& LOD, double& dpsi, double& deps, double& dx_pole, double& dy_pole, double& TAI_UTC, char interp){

    int i= -1;
    double mjd = floor(Mjd_UTC);

    for (int k = 1; k <= eop.ncols(); k++) {
        if (static_cast<int>(eop(4,k)) == static_cast<int>(mjd)) {
            i = k;
            break;
        }
    }
    if (i == -1) {
        std::cerr << "MJD not found in eop data" << std::endl;
        exit(EXIT_FAILURE);
    }

    if (interp =='l') {
        // linear interpolation
        Matrix preeop = Matrix(13,1);
        Matrix nexteop = Matrix(13,1);
        for (int j = 1; j <= 13; j++) {
            preeop(j) = eop(j,i);
            nexteop(j) = eop(j,i+1);
        }

        double mfme = 1440 * (Mjd_UTC - floor(Mjd_UTC));
        double fixf = mfme / 1440;
        // Setting of IERS Earth rotation parameters
        // (UT1-UTC [s], TAI-UTC [s], x ["], y ["])
        x_pole =    preeop(5) + (nexteop(5) - preeop(5)) * fixf;
        y_pole =    preeop(6) + (nexteop(6) - preeop(6)) * fixf;
        UT1_UTC =   preeop(7) + (nexteop(7) - preeop(7)) * fixf;
        LOD =       preeop(8) + (nexteop(8) - preeop(8)) * fixf;
        dpsi =      preeop(9) + (nexteop(9) - preeop(9)) * fixf;
        deps =      preeop(10) + (nexteop(10) - preeop(10)) * fixf;
        dx_pole =   preeop(11) + (nexteop(11) - preeop(11)) * fixf;
        dy_pole =   preeop(12) + (nexteop(12) - preeop(12)) * fixf;
        TAI_UTC =   preeop(13);

        x_pole = x_pole /Arcs;  // Pole coordinate [rad]
        y_pole = y_pole /Arcs;  // Pole coordinate [rad]
        dpsi = dpsi /Arcs;
        deps = deps /Arcs;
        dx_pole = dx_pole /Arcs; // Pole coordinate [rad]
        dy_pole = dy_pole /Arcs; // Pole coordinate [rad]
    }else if (interp =='n') {
        x_pole  = eop(5, i) / Arcs;    // Pole coordinate [rad]
        y_pole  = eop(6,i) / Arcs;    // Pole coordinate [rad]
        UT1_UTC = eop(7,i);           // UT1-UTC time difference [s]
        LOD     = eop(8,i);           // Length of day [s]
        dpsi    = eop(9,i) / Arcs;
        deps    = eop(10,i) / Arcs;
        dx_pole = eop(11,i) / Arcs;   // Pole coordinate [rad]
        dy_pole = eop(12,i) / Arcs;   // Pole coordinate [rad]
        TAI_UTC = eop(13,i);          // TAI-UTC time difference [s]
    }
}