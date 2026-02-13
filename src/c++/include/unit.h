/// @file unit.h
/// @brief Header file for unit function
/// @date 2025-04-15
/// @author Martín Hernández Tonzán

#ifndef C___UNIT_H
#define C___UNIT_H

#include "Matrix.h"

/// @ brief Calculates a unit vector given the original vector. if a zero vector is input, the vector is set to zero.
/// @param[in] vec Matrix representing vector
/// @return Matrix representing unit vector
Matrix unit(const Matrix& vec);

#endif //C___UNIT_H
