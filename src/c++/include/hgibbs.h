/// @file hgibbs.h
/// @brief Header file for hgibbs function
/// @date 2025-05-07
/// @author Martín Hernández Tonzán

#ifndef C___HGIBBS_H
#define C___HGIBBS_H

#include "Matrix.h"

/// @brief this function performs the gibbs method of orbit determination. this method determines the velocity at the middle point of the 3 given position vectors.
/// @param[in] r1 ijk position vector #1 [m]
/// @param[in] r2 ijk position vector #2 [m]
/// @param[in] r3 ijk position vector #3 [m]
/// @param[in] Mjd1 julian date of 1st sighting [days from 4713 bc]
/// @param[in] Mjd2 julian date of 2nd sighting [days from 4713 bc]
/// @param[in] Mjd3 julian date of 3rd sighting [days from 4713 bc]
/// @param[out] v2 ijk velocity vector for r2 [m/s]
/// @param[out] theta angl between vectors [rad]
/// @param[out] theta1 additional angle parameter [rad]
/// @param[out] copa additional angle parameter [rad]
/// @param[out] error flag indicating sucess 'ok', 'error', ...
void hgibbs(const Matrix& r1, const Matrix& r2, const Matrix& r3, double Mjd1, double Mjd2, double Mjd3, Matrix &v2, double &theta, double &theta1, double &copa, char*&error);

#endif //C___HGIBBS_H
