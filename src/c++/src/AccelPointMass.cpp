/// @file AccelPointMass.cpp
/// @brief Source file for AccelPointMass function
/// @date 2025-04-14
/// @author Martín Hernández Tonzán

#include "../include/AccelPointMass.h"

Matrix AccelPointMass(Matrix r, const Matrix& s, double GM){
    // Relative position vector of satellite w.r.t. point mass
    Matrix d = (r-s);

    // Acceleration
    return -GM * ( d/pow(Matrix::norm(d),3) + s/(pow(Matrix::norm(s),3)) );
}