/// @file R_y.h
/// @brief Header file for R_y function
/// @date 2025-04-09
/// @author Martín Hernández Tonzán

#ifndef C___R_Y_H
#define C___R_Y_H

#include "Matrix.h"

/// @brief Returns the 3x3 rotation matrix for a rotation around the Y-axis
/// @param[in] angle rotation angle in radians
/// @return 3x3 Matrix representing the rotation
Matrix R_y(double angle);

#endif //C___R_Y_H
