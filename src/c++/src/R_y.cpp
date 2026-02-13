/// @file R_y.cpp
/// @brief Source file for R_y function
/// @date 2025-04-09
/// @author Martín Hernández Tonzán

#include "../include/Matrix.h"
#include <cmath>

Matrix R_y(double angle){
    Matrix rotmat(3,3);
    double C = cos(angle);
    double S = sin(angle);

    rotmat(1,1) =   C;  rotmat(1,2) = 0.0;  rotmat(1,3) = -1.0*S;
    rotmat(2,1) = 0.0;  rotmat(2,2) = 1.0;  rotmat(2,3) =    0.0;
    rotmat(3,1) =   S;  rotmat(3,2) = 0.0;  rotmat(3,3) =      C;

    return rotmat;
}