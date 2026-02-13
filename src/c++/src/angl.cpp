/// @file angl.cpp
/// @brief Source file for angl function
/// @date 2025-05-07
/// @author Martín Hernández Tonzán

#include <cmath>
#include "../include/angl.h"

double sign(double temp) {
    if(temp < 0){
        return -1.0;
    } else if (temp > 0){
        return 1.0;
    } else{
        return 0.0;
    }
}

double angl (const Matrix& vec1, const Matrix& vec2){
    double small     = 0.00000001;
    double undefined = 999999.1;

    double magv1 = Matrix::norm(vec1);
    double magv2 = Matrix::norm(vec2);

    if (magv1*magv2 > pow(small,2)){
        double temp= Matrix::dot(vec1,vec2) / (magv1*magv2);
        if (fabs( temp ) > 1.0){
            temp = sign(temp) * 1.0;
        }
        return acos( temp );
    } else{
        return undefined;
    }
}