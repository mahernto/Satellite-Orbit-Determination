/// @file JPL_Eph_DE430.cpp
/// @brief Source file for JPL_Eph_DE430 function
/// @date 2025-04-16
/// @author Martín Hernández Tonzán

#include <iostream>
#include <iomanip>
#include "../include/JPL_Eph_DE430.h"
#include "../EKF_Global.h"
#include "../include/Cheb3D.h"

void JPL_Eph_DE430(double Mjd_TDB,
                   Matrix& r_Mercury, Matrix& r_Venus, Matrix& r_Earth, Matrix& r_Mars,
                   Matrix& r_Jupiter, Matrix& r_Saturn, Matrix& r_Uranus, Matrix& r_Neptune,
                   Matrix& r_Pluto, Matrix& r_Moon, Matrix& r_Sun){
    double JD = Mjd_TDB + 2400000.5;
    int fila = -1;
    for (int i = 1; i <= PC.nfils(); ++i) {
        if (PC(i, 1) <= JD && JD <= PC(i,2)) {
            fila = i;
            break;
        }
    }

    if (fila == -1) {
        std::cerr << "Date out of range in JPL_Eph_ED430.cpp" << std::endl;
        exit(EXIT_FAILURE);
    }

    Matrix PCtemp(PC.ncols(), 1);

    for (int j = 1; j <= PC.ncols(); ++j) {
        PCtemp(j) = PC(fila, j);
    }

    double t1 = PCtemp(1)-2400000.5; // MJD at start of interval

    double dt = Mjd_TDB - t1;

    Matrix temp = Matrix::range(231,270,13);
    Matrix Cx_Earth = PCtemp.subvector(temp(1), temp(2) - 1);
    Matrix Cy_Earth = PCtemp.subvector(temp(2), temp(3) - 1);
    Matrix Cz_Earth = PCtemp.subvector(temp(3), temp(4) - 1);
    temp = temp+39;
    Matrix Cx = PCtemp.subvector(temp(1), temp(2) - 1);
    Matrix Cy = PCtemp.subvector(temp(2), temp(3) - 1);
    Matrix Cz = PCtemp.subvector(temp(3), temp(4) - 1);
    Cx_Earth = Cx_Earth.concatRow(Cx);
    Cy_Earth = Cy_Earth.concatRow(Cy);
    Cz_Earth = Cz_Earth.concatRow(Cz);
    int j;
    double Mjd0;
    if (0<=dt && dt<=16){
        j=0;
        Mjd0 = t1;
    } else if(16<dt && dt<=32){
        j=1;
        Mjd0 = t1+16*j;
    }
    int index = 13 * j;
    Matrix coeffX = Cx_Earth.subvector(index + 1, index + 13);
    Matrix coeffY = Cy_Earth.subvector(index + 1, index + 13);
    Matrix coeffZ = Cz_Earth.subvector(index + 1, index + 13);
    r_Earth = 1000.0 * Cheb3D(Mjd_TDB, 13, Mjd0, Mjd0 + 16, coeffX, coeffY, coeffZ);

    temp = Matrix::range(441,480,13);
    Matrix Cx_Moon = PCtemp.subvector(temp(1), temp(2) - 1);
    Matrix Cy_Moon = PCtemp.subvector(temp(2), temp(3) - 1);
    Matrix Cz_Moon = PCtemp.subvector(temp(3), temp(4) - 1);

    for (int i = 1; i <= 7; i++){
        temp = temp+39;
        Cx = PCtemp.subvector(temp(1), temp(2) - 1);
        Cy = PCtemp.subvector(temp(2), temp(3) - 1);
        Cz = PCtemp.subvector(temp(3), temp(4) - 1);
        Cx_Moon = Cx_Moon.concatRow(Cx);
        Cy_Moon = Cy_Moon.concatRow(Cy);
        Cz_Moon = Cz_Moon.concatRow(Cz);
    }
    if (0<=dt && dt<=4){
        j=0;
        Mjd0 = t1;
    } else if(4<dt && dt<=8){
        j=1;
        Mjd0 = t1+4*j;
    } else if(8<dt && dt<=12){
        j=2;
        Mjd0 = t1+4*j;
    } else if(12<dt && dt<=16){
        j=3;
        Mjd0 = t1+4*j;
    } else if(16<dt && dt<=20){
        j=4;
        Mjd0 = t1+4*j;
    } else if(20<dt && dt<=24){
        j=5;
        Mjd0 = t1+4*j;
    } else if(24<dt && dt<=28){
        j=6;
        Mjd0 = t1+4*j;
    } else if(28<dt && dt<=32){
        j=7;
        Mjd0 = t1+4*j;
    }
    index = 13 * j;
    coeffX = Cx_Moon.subvector(index + 1, index + 13);
    coeffY = Cy_Moon.subvector(index + 1, index + 13);
    coeffZ = Cz_Moon.subvector(index + 1, index + 13);

    r_Moon = 1e3*Cheb3D(Mjd_TDB, 13, Mjd0, Mjd0+4, coeffX, coeffY, coeffZ);

    temp = Matrix::range(753,786,11);
    Matrix Cx_Sun = PCtemp.subvector(temp(1), temp(2) - 1);
    Matrix Cy_Sun = PCtemp.subvector(temp(2), temp(3) - 1);
    Matrix Cz_Sun = PCtemp.subvector(temp(3), temp(4) - 1);
    temp = temp+33;
    Cx = PCtemp.subvector(temp(1), temp(2) - 1);
    Cy = PCtemp.subvector(temp(2), temp(3) - 1);
    Cz = PCtemp.subvector(temp(3), temp(4) - 1);
    Cx_Sun = Cx_Sun.concatRow(Cx);
    Cy_Sun = Cy_Sun.concatRow(Cy);
    Cz_Sun = Cz_Sun.concatRow(Cz);
    if (0<=dt && dt<=16){
        j=0;
        Mjd0 = t1;
    } else if(16<dt && dt<=32){
        j=1;
        Mjd0 = t1+16*j;
    }
    index = 11 * j;
    coeffX = Cx_Sun.subvector(index + 1, index + 11);
    coeffY = Cy_Sun.subvector(index + 1, index + 11);
    coeffZ = Cz_Sun.subvector(index + 1, index + 11);

    r_Sun = 1e3*Cheb3D(Mjd_TDB, 11, Mjd0, Mjd0+16, coeffX, coeffY, coeffZ);

    temp = Matrix::range(3,45,14);
    Matrix Cx_Mercury = PCtemp.subvector(temp(1), temp(2) - 1);
    Matrix Cy_Mercury = PCtemp.subvector(temp(2), temp(3) - 1);
    Matrix Cz_Mercury = PCtemp.subvector(temp(3), temp(4) - 1);
    for (int i = 1; i <= 3; i++){
        temp = temp+42;
        Cx = PCtemp.subvector(temp(1), temp(2) - 1);
        Cy = PCtemp.subvector(temp(2), temp(3) - 1);
        Cz = PCtemp.subvector(temp(3), temp(4) - 1);
        Cx_Mercury = Cx_Mercury.concatRow(Cx);
        Cy_Mercury = Cy_Mercury.concatRow(Cy);
        Cz_Mercury = Cz_Mercury.concatRow(Cz);
    }

    if (0<=dt && dt<=8){
        j=0;
        Mjd0 = t1;
    } else if(8<dt && dt<=16){
        j=1;
        Mjd0 = t1+8*j;
    } else if (16<dt && dt<=24){
        j=2;
        Mjd0 = t1+8*j;
    } else if(24<dt && dt<=32){
        j=3;
        Mjd0 = t1+8*j;
    }
    index = 14 * j;
    coeffX = Cx_Mercury.subvector(index + 1, index + 14);
    coeffY = Cy_Mercury.subvector(index + 1, index + 14);
    coeffZ = Cz_Mercury.subvector(index + 1, index + 14);
    r_Mercury = 1e3*Cheb3D(Mjd_TDB, 14, Mjd0, Mjd0+8, coeffX, coeffY, coeffZ);

    temp = Matrix::range(171,201,10);
    Matrix Cx_Venus = PCtemp.subvector(temp(1), temp(2) - 1);
    Matrix Cy_Venus = PCtemp.subvector(temp(2), temp(3) - 1);
    Matrix Cz_Venus = PCtemp.subvector(temp(3), temp(4) - 1);
    temp = temp+30;
    Cx = PCtemp.subvector(temp(1), temp(2) - 1);
    Cy = PCtemp.subvector(temp(2), temp(3) - 1);
    Cz = PCtemp.subvector(temp(3), temp(4) - 1);
    Cx_Venus = Cx_Venus.concatRow(Cx);
    Cy_Venus = Cy_Venus.concatRow(Cy);
    Cz_Venus = Cz_Venus.concatRow(Cz);
    if (0<=dt && dt<=16){
        j=0;
        Mjd0 = t1;
    } else if(16<dt && dt<=32){
        j=1;
        Mjd0 = t1+16*j;
    }
    index = 10 * j;
    coeffX = Cx_Venus.subvector(index + 1, index + 10);
    coeffY = Cy_Venus.subvector(index + 1, index + 10);
    coeffZ = Cz_Venus.subvector(index + 1, index + 10);
    r_Venus = 1e3*Cheb3D(Mjd_TDB, 10, Mjd0, Mjd0+16, coeffX, coeffY, coeffZ);

    temp = Matrix::range(309,342,11);
    Matrix Cx_Mars = PCtemp.subvector(temp(1), temp(2) - 1);
    Matrix Cy_Mars = PCtemp.subvector(temp(2), temp(3) - 1);
    Matrix Cz_Mars = PCtemp.subvector(temp(3), temp(4) - 1);
    j=0;
    Mjd0 = t1;
    index = 11 * j;
    coeffX = Cx_Mars.subvector(index + 1, index + 11);
    coeffY = Cy_Mars.subvector(index + 1, index + 11);
    coeffZ = Cz_Mars.subvector(index + 1, index + 11);
    r_Mars = 1e3*Cheb3D(Mjd_TDB, 11, Mjd0, Mjd0+32, coeffX, coeffY, coeffZ);

    temp = Matrix::range(342,366,8);
    Matrix Cx_Jupiter = PCtemp.subvector(temp(1), temp(2) - 1);
    Matrix Cy_Jupiter = PCtemp.subvector(temp(2), temp(3) - 1);
    Matrix Cz_Jupiter = PCtemp.subvector(temp(3), temp(4) - 1);
    j=0;
    Mjd0 = t1;
    index = 8 * j;
    coeffX = Cx_Jupiter.subvector(index + 1, index + 8);
    coeffY = Cy_Jupiter.subvector(index + 1, index + 8);
    coeffZ = Cz_Jupiter.subvector(index + 1, index + 8);
    r_Jupiter = 1e3*Cheb3D(Mjd_TDB, 8, Mjd0, Mjd0+32, coeffX, coeffY, coeffZ);

    temp = Matrix::range(366,387,7);
    Matrix Cx_Saturn = PCtemp.subvector(temp(1), temp(2) - 1);
    Matrix Cy_Saturn = PCtemp.subvector(temp(2), temp(3) - 1);
    Matrix Cz_Saturn = PCtemp.subvector(temp(3), temp(4) - 1);
    j=0;
    Mjd0 = t1;
    index = 7 * j;
    coeffX = Cx_Saturn.subvector(index + 1, index + 7);
    coeffY = Cy_Saturn.subvector(index + 1, index + 7);
    coeffZ = Cz_Saturn.subvector(index + 1, index + 7);
    r_Saturn = 1e3*Cheb3D(Mjd_TDB, 7, Mjd0, Mjd0+32, coeffX, coeffY, coeffZ);

    temp = Matrix::range(387,405,6);
    Matrix Cx_Uranus = PCtemp.subvector(temp(1), temp(2) - 1);
    Matrix Cy_Uranus = PCtemp.subvector(temp(2), temp(3) - 1);
    Matrix Cz_Uranus = PCtemp.subvector(temp(3), temp(4) - 1);
    j=0;
    Mjd0 = t1;
    index = 6 * j;
    coeffX = Cx_Uranus.subvector(index + 1, index + 6);
    coeffY = Cy_Uranus.subvector(index + 1, index + 6);
    coeffZ = Cz_Uranus.subvector(index + 1, index + 6);
    r_Uranus = 1e3*Cheb3D(Mjd_TDB, 6, Mjd0, Mjd0+32, coeffX, coeffY, coeffZ);

    temp = Matrix::range(405,423,6);
    Matrix Cx_Neptune = PCtemp.subvector(temp(1), temp(2) - 1);
    Matrix Cy_Neptune = PCtemp.subvector(temp(2), temp(3) - 1);
    Matrix Cz_Neptune = PCtemp.subvector(temp(3), temp(4) - 1);
    j=0;
    Mjd0 = t1;
    index = 6 * j;
    coeffX = Cx_Neptune.subvector(index + 1, index + 6);
    coeffY = Cy_Neptune.subvector(index + 1, index + 6);
    coeffZ = Cz_Neptune.subvector(index + 1, index + 6);
    r_Neptune = 1e3*Cheb3D(Mjd_TDB, 6, Mjd0, Mjd0+32, coeffX, coeffY, coeffZ);

    temp = Matrix::range(423,441,6);
    Matrix Cx_Pluto = PCtemp.subvector(temp(1), temp(2) - 1);
    Matrix Cy_Pluto = PCtemp.subvector(temp(2), temp(3) - 1);
    Matrix Cz_Pluto = PCtemp.subvector(temp(3), temp(4) - 1);
    j=0;
    Mjd0 = t1;
    index = 6 * j;
    coeffX = Cx_Pluto.subvector(index + 1, index + 6);
    coeffY = Cy_Pluto.subvector(index + 1, index + 6);
    coeffZ = Cz_Pluto.subvector(index + 1, index + 6);
    r_Pluto = 1e3*Cheb3D(Mjd_TDB, 6, Mjd0, Mjd0+32, coeffX, coeffY, coeffZ);

    temp = Matrix::range(819,839,10);
    Matrix Cx_Nutations = PCtemp.subvector(temp(1), temp(2) - 1);
    Matrix Cy_Nutations = PCtemp.subvector(temp(2), temp(3) - 1);
    for (int i=1; i <=3; i++){
        temp = temp+20;
        Cx = PCtemp.subvector(temp(1), temp(2) - 1);
        Cy = PCtemp.subvector(temp(2), temp(3) - 1);
        Cx_Nutations = Cx_Nutations.concatRow(Cx);
        Cy_Nutations = Cy_Nutations.concatRow(Cy);
    }
    if (0<=dt && dt<=8){
        j=0;
        Mjd0 = t1;
    } else if(8<dt && dt<=16){
        j=1;
        Mjd0 = t1+8*j;
    } else if (16<dt && dt<=24){
        j=2;
        Mjd0 = t1+8*j;
    } else if(24<dt && dt<=32){
        j=3;
        Mjd0 = t1+8*j;
    }
    index = 10 * j;
    coeffX = Cx_Nutations.subvector(index + 1, index + 10);
    coeffY = Cy_Nutations.subvector(index + 1, index + 10);
    Matrix Nutations = Cheb3D(Mjd_TDB, 10, Mjd0, Mjd0+8, coeffX, coeffY,Matrix(10,1));

    temp = Matrix::range(899,929,10);
    Matrix Cx_Librations = PCtemp.subvector(temp(1), temp(2) - 1);
    Matrix Cy_Librations = PCtemp.subvector(temp(2), temp(3) - 1);
    Matrix Cz_Librations = PCtemp.subvector(temp(2), temp(4) - 1);
    for (int i=1; i <= 3; i++){
        temp = temp+30;
        Cx = PCtemp.subvector(temp(1), temp(2) - 1);
        Cy = PCtemp.subvector(temp(2), temp(3) - 1);
        Cz = PCtemp.subvector(temp(3), temp(4) - 1);
        Cx_Librations = Cx_Librations.concatRow(Cx);
        Cy_Librations = Cy_Librations.concatRow(Cy);
        Cz_Librations = Cz_Librations.concatRow(Cz);
    }
    if (0<=dt && dt<=8){
        j=0;
        Mjd0 = t1;
    } else if(8<dt && dt<=16){
        j=1;
        Mjd0 = t1+8*j;
    } else if (16<dt && dt<=24){
        j=2;
        Mjd0 = t1+8*j;
    } else if(24<dt && dt<=32){
        j=3;
        Mjd0 = t1+8*j;
    }
    index = 10 * j;
    coeffX = Cx_Librations.subvector(index + 1, index + 10);
    coeffY = Cy_Librations.subvector(index + 1, index + 10);
    coeffZ = Cz_Librations.subvector(index + 1, index + 10);
    Matrix Librations = Cheb3D(Mjd_TDB, 10, Mjd0, Mjd0+8, coeffX, coeffY, coeffZ);
    double EMRAT = 81.30056907419062; // DE430
    double EMRAT1 = 1/(1+EMRAT);
    r_Earth = r_Earth-EMRAT1*r_Moon;
    r_Mercury = -r_Earth+r_Mercury;
    r_Venus = -r_Earth+r_Venus;
    r_Mars = -r_Earth+r_Mars;
    r_Jupiter = -r_Earth+r_Jupiter;
    r_Saturn = -r_Earth+r_Saturn;
    r_Uranus = -r_Earth+r_Uranus;
    r_Neptune = -r_Earth+r_Neptune;
    r_Pluto = -r_Earth+r_Pluto;
    r_Sun = -r_Earth+r_Sun;
}