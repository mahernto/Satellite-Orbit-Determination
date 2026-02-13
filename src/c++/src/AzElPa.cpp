/// @file AzElPa.cpp
/// @brief Source file for AzElPa function
/// @date 2025-04-14
/// @author Martín Hernández Tonzán

#include "../include/AzElPa.h"
#include "../include/SAT_Const.h"

void AzElPa(const Matrix& s, double &Az, double &El, Matrix &dAds, Matrix &dEds){
    double rho = sqrt(s(1)*s(1)+s(2)*s(2));

    // Angles
    Az = atan2(s(1),s(2));

    if (Az<0.0)
        Az = Az+pi2;

    El = atan ( s(3) / rho );

    // Partials
    dAds(1) = s(2)/(rho*rho);
    dAds(2) = -s(1)/(rho*rho);
    dAds(3) = 0.0;
    dEds(1) = (-s(1)*s(3)/rho) / Matrix::dot(s,s);
    dEds(2) = (-s(2)*s(3)/rho) / Matrix::dot(s,s);
    dEds(3) =  rho / Matrix::dot(s,s);

}