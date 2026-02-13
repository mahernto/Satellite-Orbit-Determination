/// @file Frac.cpp
/// @brief Source file for Frac function
/// @date 2025-04-10
/// @author Martín Hernández Tonzán

#include <cmath>
#include "../include/Frac.h"

double Frac(double x){
    return x - floor(x);
}