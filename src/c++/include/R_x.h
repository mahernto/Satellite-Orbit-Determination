/// @file R_x.h
/// @brief Header file for R_x function
/// @date 2025-04-09
/// @author Martín Hernández Tonzán

#ifndef C___R_X_H
#define C___R_X_H

#include "Matrix.h"
#include <valarray>

/// @brief Returns the 3x3 rotation matrix for a rotation around the X-axis
/// @param[in] angle rotation angle in radians
/// @return 3x3 Matrix representing the rotation
Matrix R_x(double angle);


#endif //C___R_X_H
