//
// Created by martin on 08/05/2025.
//

#include <cmath>
#include "vector.h"

void cross(const double a[3], const double b[3], double result[3]) {
    result[0] = a[1]*b[2] - a[2]*b[1];
    result[1] = a[2]*b[0] - a[0]*b[2];
    result[2] = a[0]*b[1] - a[1]*b[0];
}

// Producto punto
double dot(const double a[3], const double b[3]) {
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}

// Norma de vector 3D
double norm(const double v[3]) {
    return sqrt(dot(v, v));
}