/// @file EKF_Global.h
/// @brief Header file for global variables in the project
/// @date 2025-04-09
/// @version 1.0

#ifndef EKF_TEST_H
#define EKF_TEST_H

#include "include/Matrix.h"

extern Matrix eop;
extern Matrix PC;
extern Matrix Cnm;
extern Matrix Snm;
extern Matrix obs;
struct AuxParam {
    double Mjd_UTC;
    int n;
    int m;
    bool sun;
    bool moon;
    bool planets;
    double Mjd_TT;
};
extern AuxParam auxParam;

#endif // EKF_TEST_H
