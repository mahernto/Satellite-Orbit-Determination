/// @file R_x.cpp
/// @brief Source file for R_x function
/// @date 2025-04-09
/// @author Martín Hernández Tonzán

#include "../include/R_x.h"

Matrix R_x(double angle){
    Matrix rotmat(3,3);
    double C = cos(angle);
    double S = sin(angle);

    rotmat(1,1) = 1.0;  rotmat(1,2) =    0.0;  rotmat(1,3) = 0.0;
    rotmat(2,1) = 0.0;  rotmat(2,2) =      C;  rotmat(2,3) =   S;
    rotmat(3,1) = 0.0;  rotmat(3,2) = -1.0*S;  rotmat(3,3) =   C;

    return rotmat;
}