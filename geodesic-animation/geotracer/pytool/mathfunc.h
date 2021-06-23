#ifndef MATHFUNC_H
#define MATHFUNC_H

#include <stdio.h>
#include <math.h>


int invert_44(double gcov[4][4], double gcon[4][4]);

// borrowed from HARM
int invert_matrix(double Am[][4], double Aminv[][4]);
int LU_decompose(double A[][4], int permute[]);
void LU_substitution(double A[][4], double B[], int permute[]);

#endif // MATHFUNC_H