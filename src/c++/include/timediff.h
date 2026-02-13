/// @file timediff.h
/// @brief Header file for timediff function
/// @date 2025-04-14
/// @author Martín Hernández Tonzán

#ifndef C___TIMEDIFF_H
#define C___TIMEDIFF_H

/// @brief Calculates various time differences: UT1-TAI, UTC-GPS, UT1-GPS, TT-UTC, and GPS-UTC
/// @param[in] UT1_UTC The difference between UT1 and UTC in seconds
/// @param[in] TAI_UTC The difference between TAI and UTC in seconds
/// @param[out] UT1_TAI The difference between UT1 and TAI in seconds
/// @param[out] UTC_GPS The difference between UTC and GPS time in seconds
/// @param[out] UT1_GPS The difference between UT1 and GPS time in seconds
/// @param[out] TT_UTC The difference between TT (Terrestrial Time) and UTC in seconds
/// @param[out] GPS_UTC The difference between GPS time and UTC in seconds
void timediff(double UT1_UTC, double TAI_UTC, double &UT1_TAI, double &UTC_GPS, double &UT1_GPS, double &TT_UTC, double &GPS_UTC);

#endif //C___TIMEDIFF_H
