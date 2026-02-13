/// @file TimeUpdate.h
/// @brief Header file for TimeUpdate function
/// @date 2025-04-14
/// @author Martín Hernández Tonzán

#ifndef C___TIMEUPDATE_H
#define C___TIMEUPDATE_H

#include "Matrix.h"

/// @brief Performs the time update step of the Kalman filter
/// @param[in] Phi State transition matrix
/// @param[out] P Covariance matrix of the estimation error
/// @param[in] Qdt Process noise covariance matrix (or scalar) multiplied by time step. Defaults to 0 if not provided
Matrix TimeUpdate(const Matrix& P, Matrix Phi, const Matrix& Qdt= Matrix(6,6));

#endif //C___TIMEUPDATE_H
