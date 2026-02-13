/// @file AzElPa.h
/// @brief Header file for AzElPa function
/// @date 2025-04-14
/// @author Martín Hernández Tonzán

#ifndef C___AZELPA_H
#define C___AZELPA_H

#include "Matrix.h"

/// @ brief Computes azimuth, elevation and partials from local tangent coordinates
/// @param[in] s Topocentric local tangent coordinates (East-North-Zenith frame)
/// @param[out] Az Azimuth [rad]
/// @param[out] El Elevation [rad]
/// @param[out] dAds Partials of azimuth w.r.t. s
/// @param[out] dEds Partials of elevation w.r.t. s
void AzElPa(const Matrix& s, double &Az, double &El, Matrix &dAds, Matrix &dEds);

#endif //C___AZELPA_H
