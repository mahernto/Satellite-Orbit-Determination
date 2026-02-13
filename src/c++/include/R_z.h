/// @file R_z.h
/// @brief Header file for R_z function
/// @date 2025-04-09
/// @author Martín Hernández Tonzán

#ifndef C___R_Z_H
#define C___R_Z_H

#include "Matrix.h"

/// @brief Returns the 3x3 rotation matrix for a rotation around the Z-axis
/// @param[in] angle rotation angle in radians
/// @return 3x3 Matrix representing the rotation
Matrix R_z(double angle);

#endif //C___R_Z_H
