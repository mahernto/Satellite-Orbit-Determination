/// @file MeasUpdate.h
/// @brief Header file for MeasUpdate function
/// @date 2025-04-14
/// @author Martín Hernández Tonzán

#ifndef C___MEASUPDATE_H
#define C___MEASUPDATE_H

#include "Matrix.h"

/// @brief Performs the measurement update step of the Kalman Filter.
/// @param[in] x Prior state estimate vector
/// @param[in] z Measurement vector
/// @param[in] g Expected measurement (predicted from state)
/// @param[in] s Standard deviations of measurement errors (sqrt of diagonal of R)
/// @param[in] G Measurement matrix (Jacobian of g w.r.t. state)
/// @param[in] P Prior covariance matrix
/// @param[in] n Dimension of the state vector
/// @param[in] n Dimension of the measurement vector
/// @param[out] K Kalman gain matrix
/// @param[out] Y Updated state estimate vector
/// @param[out] Pout Updated covariance matrix
void MeasUpdate(Matrix x, Matrix z, const Matrix& g, const Matrix& s, Matrix G, Matrix P, int n, Matrix &K, Matrix &Y, Matrix &Pout);

#endif //C___MEASUPDATE_H
