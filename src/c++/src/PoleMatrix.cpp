/// @file PoleMatrix.cpp
/// @brief Source file for PoleMatrix function
/// @date 2025-04-16
/// @author Martín Hernández Tonzán

#include "../include/PoleMatrix.h"
#include "../include/R_x.h"
#include "../include/R_y.h"

Matrix PoleMatrix(double xp, double yp){
    return R_y(-xp) * R_x(-yp);
}