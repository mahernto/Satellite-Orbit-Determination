/// @file GHAMatrix.cpp
/// @brief Source file for GHAMatrix function
/// @date 2025-04-16
/// @author Martín Hernández Tonzán

#include "../include/GHAMatrix.h"
#include "../include/gast.h"
#include "../include/R_z.h"

Matrix GHAMatrix(double Mjd_UT1){
    return R_z( gast(Mjd_UT1) );
}