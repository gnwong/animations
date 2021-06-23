#ifndef GEOMETRY_H
#define GEOMETRY_H

#include <math.h>

// if (1), central pixel will have k_phi = 0, 
// otherwise will use old method
#define USE_IMPACT_ZEROING (1)

// set quality factor
void set_eps(double eps);

void get_connection(double X[4], double lconn[4][4][4]);

double stepsize(double X[4], double Kcon[4]);
void push_photon(double X[4], double Kcon[4], double dl, double Xhalf[4], double Kconhalf[4]);

double gdet_func(double gcov[][4]);

// tetrad business
void tetrad_to_coordinate(double Econ[4][4], double K_tetrad[4], double K[4]);
void set_Econ_from_trial(double Econ[4], int defdir, double trial[4]);
void project_out(double vcona[4], double vconb[4], double Gcov[4][4]);
void normalize(double vcon[4], double Gcov[4][4]);
void lower(double ucon[4], double Gcov[4][4], double ucov[4]);
void make_plasma_tetrad(double Ucon[4], double Kcon[4],
      double Bcon[4], double Gcov[4][4],
      double Econ[4][4], double Ecov[4][4]);
void make_camera_tetrad(double X[4], double Econ[4][4], double Ecov[4][4]);
double check_handedness(double Econ[4][4], double Gcov[4][4]);

#endif // GEOMETRY_H
