/// @file Mjday.h
/// @brief Header file for Mjday function
/// @date 2025-04-09
/// @author Martín Hernández Tonzán

#ifndef C___MJDAY_H
#define C___MJDAY_H

/// @brief Computes Modified Julian Date from calendar date and time
/// @param[in] yr  year
/// @param[in] mon month
/// @param[in] day day of the month
/// @param[in] hr  hour of the day [default = 0.0]
/// @param[in] min minute [default = 0.0]
/// @param[in] sec second [default = 0.0]
/// @return Modified Julian Date as a double
double Mjday(double yr, double mon, double day, double hr=0.0, double min=0.0, double sec=0.0);


#endif //C___MJDAY_H
