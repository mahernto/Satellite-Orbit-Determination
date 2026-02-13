    /// @file MeasUpdate.cpp
    /// @brief Source file for MeasUpdate function
    /// @date 2025-04-14
    /// @author Martín Hernández Tonzán

    #include <iostream>
    #include "../include/MeasUpdate.h"

    void MeasUpdate(Matrix x, Matrix z, const Matrix& g, const Matrix& s, Matrix G, Matrix P, int n, Matrix &K, Matrix &Y, Matrix &Pout){
        int m = z.nfils();

        Matrix Inv_W(m,m);

        for (int i=1; i <= m; i++){
            Inv_W(i,i) = s(i)*s(i); // Inverse weight (measurement covariance)
        }
        // Kalman gain
        K = P * G.transpose() * (Inv_W + G * P * G.transpose()).inv();

        // State update
        Y = x + (K*(z-g));

        // Covariance update
        Pout = (Matrix::eye(n)-K*G)*P;
    }