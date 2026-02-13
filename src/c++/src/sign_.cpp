/// @file sign_.cpp
/// @brief Source file for sign_ function
/// @date 2025-04-16
/// @author Martín Hernández Tonzán

#include "../include/sign_.h"
#include <cmath>

double sign_(double a, double b){
    if (b>=0.0)
        return fabs(a);
    else
        return - fabs(a);
}