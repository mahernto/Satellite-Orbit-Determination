/// @file Mjday_TDB.h
/// @brief Header file for Mjday_TDB function
/// @date 2025-04-14
/// @author Martín Hernández Tonzán

#ifndef C___MJDAY_TDB_H
#define C___MJDAY_TDB_H

#include "Mjday.h"

/// @ brief Computes the Modified Julian Date for barycentric dynamical time
/// @param[in] Mjd_TT Modified julian date (Terrestrial Time)
/// @return Modified julian date (TDB)
double Mjday_TDB(double Mjd_TT);

#endif //C___MJDAY_TDB_H
