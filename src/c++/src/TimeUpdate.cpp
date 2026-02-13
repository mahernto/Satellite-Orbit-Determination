/// @file TimeUpdate.cpp
/// @brief Source file for TimeUpdate function
/// @date 2025-04-14
/// @author Martín Hernández Tonzán

#include "../include/TimeUpdate.h"

Matrix TimeUpdate(const Matrix& P, Matrix Phi, const Matrix& Qdt){
    return Phi*P*Phi.transpose() + Qdt;
}
