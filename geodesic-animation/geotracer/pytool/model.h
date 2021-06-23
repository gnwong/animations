#ifndef MODEL_H
#define MODEL_H

#include <stdlib.h>
#include <math.h>

int valid_geodesic(double X[4], double Kcon[4]);
void set_gcov(double X[4], double gcov[4][4]);
double root_find(double x[4]);
double theta_rootfind(double x[4]);
void set_params(char *argv[]);
void init_model(double bhspin, double input_rmax_geo);

void gcov_func(double X[4], double gcov[4][4]);
void gcov_ks(double r, double th, double gcov[4][4]);
void bl_coord(double X[4], double *r, double *th);
void bl_coord_vec(double X[4], double Xbl[4]);
void bl_coord_vec_many(double X[][4], double Xbl[][4], int n);
void set_dxdX(double X[4], double dxdX[4][4]);


#endif // MODEL_H
