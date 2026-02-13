/// @file VarEqn.h
/// @brief Header file for VarEqn function
/// @date 2025-04-16
/// @author Martín Hernández Tonzán

#ifndef C___VAREQN_H
#define C___VAREQN_H

#include "Matrix.h"

/// @brief Computes the variational equations, i.e. the derivative of the state vector and the state transition matrix
/// @param[in] x Time since epoch in [s]
/// @param[in] yPhi (6+36)-dim vector comprising the state vector (y) and the state transition matrix (Phi) in column wise storage order
Matrix VarEqn(double x, const Matrix& yPhi);

#endif //C___VAREQN_H
