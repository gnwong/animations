
#include "model.h"

double poly_alpha = 14.;
double mks_smooth = 0.5;
double hslope = 0.3;
double poly_xt = 0.82;
double startx[4] = { 0., 0.19685238361579208, 0, 0, }; // confirmed Ma+0.94 384 2020.06.27 gnw

double a = 0.;

double rmax_geo = 100.;



double theta_func(double X[4])
{
  double r, th;
  bl_coord(X, &r, &th);
  return th;
}

double theta_rootfind(double x[4]) 
{
  return root_find(x);
}

double root_find(double x[4])
{
  double SMALL = 1.e-40;


    double th = x[2];
    double thb, thc;
    double theta_func(double X[4]);

    double Xa[4], Xb[4], Xc[4];
    Xa[1] = log(x[1]);
    Xa[3] = x[3];
    Xb[1] = Xa[1];
    Xb[3] = Xa[3];
    Xc[1] = Xa[1];
    Xc[3] = Xa[3];

    if (x[2] < M_PI / 2.) {
      Xa[2] = 0. - SMALL;
      Xb[2] = 0.5 + SMALL;
    } else {
      Xa[2] = 0.5 - SMALL;
      Xb[2] = 1. + SMALL;
    }

    //tha = theta_func(Xa);
    thb = theta_func(Xb);

    /* bisect for a bit */
    double tol = 1.e-6;
    for (int i = 0; i < 100; i++) {
      Xc[2] = 0.5 * (Xa[2] + Xb[2]);
      thc = theta_func(Xc);

      if ((thc - th) * (thb - th) < 0.)
        Xa[2] = Xc[2];
      else
        Xb[2] = Xc[2];

      double err = theta_func(Xc) - th;
      if (fabs(err) < tol) break;
    }

    return (Xa[2]);
}

void set_gcov(double X[4], double gcov[4][4]) {

  for (int i=0; i<4; ++i) {
    for (int j=0; j<4; ++j) {
      gcov[i][j] = 0.;
      if (i==0 && j==0) gcov[i][j] = -1.;
      else if (i == j) gcov[i][j] = 1.;
    }
  }

  gcov_func(X, gcov);
}

#include <stdio.h>
int valid_geodesic(double X[4], double Kcon[4]) {

  double reh = 1. + sqrt(1. - a*a);

  if (X[1] > log(1.1*rmax_geo) && Kcon[1] < 0.) return 0; // out too far
  if (X[1] < log(1.05*reh)) return 0.; // we've fallen into the hole!

  return 1;
}

void set_params(char *argv[]) {
  a = atof(argv[2]);
}

void init_model(double bhspin, double input_rmax_geo) {
  a = bhspin;
  rmax_geo = input_rmax_geo;
}

void gcov_ks(double r, double th, double gcov[4][4])
{
  double cth = cos(th);
  double sth = sin(th);

  double s2 = sth*sth;
  double rho2 = r*r + a*a*cth*cth;
 
  for (int mu=0; mu<4; ++mu) {
    for (int nu=0; nu<4; ++nu) {
      gcov[mu][nu] = 0.; 
    }
  }

  // compute ks metric for ks coordinates (cyclic in t,phi)
  gcov[0][0] = -1. + 2.*r/rho2;
  gcov[0][1] = 2.*r/rho2;
  gcov[0][3] = -2.*a*r*s2/rho2;

  gcov[1][0] = gcov[0][1];
  gcov[1][1] = 1. + 2.*r/rho2;
  gcov[1][3] = -a*s2*(1. + 2.*r/rho2);

  gcov[2][2] = rho2;

  gcov[3][0] = gcov[0][3];
  gcov[3][1] = gcov[1][3];
  gcov[3][3] = s2*(rho2 + a*a*s2*(1. + 2.*r/rho2));
}

void bl_coord_vec_many(double X[][4], double Xbl[][4], int n)
{
  for (int i=0; i<n; ++i) {
    bl_coord_vec(X[i], Xbl[i]);
  }
}

void bl_coord_vec(double X[4], double Xbl[4])
{
  double r, th;
  bl_coord(X, &r, &th);
  Xbl[0] = X[0];
  Xbl[1] = r;
  Xbl[2] = th;
  Xbl[3] = X[3];
}

void bl_coord(double X[4], double *r, double *th)
{
  *r = exp(X[1]);

  double poly_norm = 0.5*M_PI*1./(1. + 1./(poly_alpha + 1.)*1./pow(poly_xt, poly_alpha));


  double thG = M_PI*X[2] + ((1. - hslope)/2.)*sin(2.*M_PI*X[2]);
    double y = 2*X[2] - 1.;
    double thJ = poly_norm*y*(1. + pow(y/poly_xt,poly_alpha)/(poly_alpha+1.)) + 0.5*M_PI;
    *th = thG + exp(mks_smooth*(startx[1] - X[1]))*(thJ - thG);

}

void set_dxdX(double X[4], double dxdX[4][4])
{
  // Jacobian with respect to KS basis where X is given in
  // non-KS basis
  for (int mu=0; mu<4; ++mu) {
    for (int nu=0; nu<4; ++nu) {
      dxdX[mu][nu] = 0.;
    }
  }

  double poly_norm = 0.5*M_PI*1./(1. + 1./(poly_alpha + 1.)*1./pow(poly_xt, poly_alpha));

  // mmks
  dxdX[0][0] = 1.;
  dxdX[1][1] = exp(X[1]);
  dxdX[2][1] = -exp(mks_smooth*(startx[1]-X[1]))*mks_smooth*(
    M_PI/2. -
    M_PI*X[2] +
    poly_norm*(2.*X[2]-1.)*(1+(pow((-1.+2*X[2])/poly_xt,poly_alpha))/(1 + poly_alpha)) -
    1./2.*(1. - hslope)*sin(2.*M_PI*X[2])
    );
  dxdX[2][2] = M_PI + (1. - hslope)*M_PI*cos(2.*M_PI*X[2]) +
    exp(mks_smooth*(startx[1]-X[1]))*(
      -M_PI +
      2.*poly_norm*(1. + pow((2.*X[2]-1.)/poly_xt,poly_alpha)/(poly_alpha+1.)) +
      (2.*poly_alpha*poly_norm*(2.*X[2]-1.)*pow((2.*X[2]-1.)/poly_xt,poly_alpha-1.))/((1.+poly_alpha)*poly_xt) -
      (1.-hslope)*M_PI*cos(2.*M_PI*X[2])
      );
  dxdX[3][3] = 1.;
}

void gcov_func(double X[4], double gcov[4][4])
{
  // returns g_{munu} at location specified by X

  // despite the name, get equivalent values for
  // r, th for kerr coordinate system
  double r, th;
  bl_coord(X, &r, &th);

  // compute ks metric
  gcov_ks(r, th, gcov);

  // convert from ks metric to mks/mmks
  double dxdX[4][4];
  set_dxdX(X, dxdX);

  double gcov_ks[4][4];
  
  for (int mu=0; mu<4; ++mu) {
    for (int nu=0; nu<4; ++nu) {
      gcov_ks[mu][nu] = gcov[mu][nu];
      gcov[mu][nu] = 0.;
    }
  }

  for (int mu=0; mu<4; ++mu) {
    for (int nu=0; nu<4; ++nu) {
      for (int lam=0; lam<4; ++lam) {
        for (int kap=0; kap<4; ++kap) {
          gcov[mu][nu] += gcov_ks[lam][kap]*dxdX[lam][mu]*dxdX[kap][nu];
        }
      }
    }
  }
}

