/// @file gibbs.cpp
/// @brief Source file for gibbs function
/// @date 2025-05-07
/// @author Martín Hernández Tonzán

#include <cmath>
#include "../include/gibbs.h"
#include "../include/unit.h"
#include "../include/angl.h"
#include "../include/SAT_Const.h"

void gibbs(const Matrix& r1, const Matrix& r2, const Matrix& r3, Matrix &v2, double &theta, double &theta1, double &copa, char*&error) {
    double small = 0.00000001;
    theta = 0.0;
    error = "ok";
    theta1 = 0.0;

    double magr1 = Matrix::norm(r1);
    double magr2 = Matrix::norm(r2);
    double magr3 = Matrix::norm(r3);
    for (int i = 1; i <= 3; i++) {
        v2(i) = 0.0;
    }

    Matrix p = Matrix::cross(r2, r3);
    Matrix q = Matrix::cross(r3, r1);
    Matrix w = Matrix::cross(r1, r2);
    Matrix pn = unit(p);
    Matrix r1n = unit(r1);
    copa = asin(Matrix::dot(pn, r1n));

    if (fabs(Matrix::dot(r1n, pn)) > 0.017452406){
        error = "not coplanar";
    }

    Matrix d = p + q + w;
    double magd = Matrix::norm(d);
    Matrix n = magr1*p + magr2*q + magr3*w;
    double magn = Matrix::norm(n);
    Matrix nn = unit( n );
    Matrix dn = unit( d );

    // ------------------------------------------------------------- determine if  the orbit is possible. both d and n must be in the same direction, and non-zero. -------------------------------------------------------------
    if ( ( fabs(magd)<small ) || ( fabs(magn)<small ) || ( Matrix::dot(nn,dn) < small ) ){
        error= "impossible";
    } else{
        theta  = angl( r1,r2 );
        theta1 = angl( r2,r3 );

        // ----------- perform gibbs method to find v2 -----------
        double r1mr2= magr1-magr2;
        double r3mr1= magr3-magr1;
        double r2mr3= magr2-magr3;
        Matrix s  = r1mr2*r3 + r3mr1*r2 + r2mr3*r1;
        Matrix b  = Matrix::cross( d,r2 );
        double l  = sqrt(GM_Earth / (magd*magn) );
        double tover2= l / magr2;
        v2 = tover2 * b + l * s;
    }
}