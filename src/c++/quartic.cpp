#include <math.h>
//#include <mathbits/quartic.h>
#include "./include/rpoly.h"


#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

int main()
{
    double coef[10], zeror[10], zeroi[10];
    int i, nr, degree;
    
    coef[0] = 1.0;
    coef[1] = 0.0;
    coef[2] = -73120740632072.03125;
    coef[3] = 0.0;
    coef[4] = 0.0;
    coef[5] = -1587936795676375129762729592392515584.0;
    coef[6] = 0.0;
    coef[7] = 0.0;
    coef[8] = -11985384853690913913633237309023581769694874388537465110528.0;
    degree =8;
    
    nr = real_poly_roots(coef, degree, zeror, zeroi);

    for (i=0; i<nr; ++i) {
        printf("Raiz %d   preal = %.20lf  pimag =  %.20lf\n", i+1, zeror[i], zeroi[i]);
    }
    printf("Raices reales\n");
    for (i=0; i<nr; ++i) {
        if(fabs(zeroi[i]) < 10e-12)
            printf("Raiz %d   preal = %.20lf  pimag =  %.20lf\n", i+1, zeror[i], zeroi[i]);
    }

    return 0;
}


// FINDS THE ZEROS OF A REAL POLYNOMIAL
// OP  - DOUBLE PRECISION VECTOR OF COEFFICIENTS IN
//       ORDER OF DECREASING POWERS.
// DEGREE   - INTEGER DEGREE OF POLYNOMIAL.
// ZEROR, ZEROI - OUTPUT DOUBLE PRECISION VECTORS OF
//                REAL AND IMAGINARY PARTS OF THE
//                ZEROS.
// RETURNS NUMBER OF ZEROS FOUND.
