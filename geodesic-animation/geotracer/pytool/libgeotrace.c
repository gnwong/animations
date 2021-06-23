#include <stdio.h>
#include <string.h>

#include "mathfunc.h"
#include "geometry.h"
#include "model.h"

#define MAX_GEODESIC_STEPS (1000000)


void calc (const double *indata, size_t size, double *outdata) {
  for (size_t i=0; i<size; ++i) {
    outdata[i] = indata[i] * 2.0;
  }
}


/* construct orthonormal tetrad.
 *   e^0 along Ucam
 *   e^3 outward (!) along radius vector
 *   e^2 toward north pole of coordinate system
 *   ("y" for the image plane)
 *   e^1 in the remaining direction
 *   ("x" for the image plane)
 */
void init_XK (int i, int j, double di, double dj, double Xcam[4], double fovx, double fovy, double X[4], double Kcon[4], double Kcov[4], int NX, int NY) {

  #define P4VECHP(S,X) fprintf(stderr, "%s: %.17g %.17g %.17g %.17g\n", S, X[0], X[1], X[2], X[3]);

  double Econ[4][4];
  double Ecov[4][4];
  double Kcon_tetrad[4];
  
  {
    // make_camera_tetrad copied directly here. this is CAMERA_CENTER_ZAMO  
    double Gcov[4][4], Gcon[4][4];
    gcov_func(Xcam, Gcov);
    invert_44(Gcov, Gcon);
    double Ucam[4], Kcon[4], trial[4];

    // center the camera according to impact parameter, i.e., make it
    // so that Kcontetrad = ( 1, 0, 0, 1 ) corresponds to an outgoing
    // wavevector with zero angular momentum / zero impact parameter.

    // use normal observer velocity. this forces (Gcov.Econ[0])[3] = 0.
    trial[0] = -1.;
    trial[1] = 0.;
    trial[2] = 0.;
    trial[3] = 0.;
    lower(trial, Gcon, Ucam);

    // set Kcon (becomes Econ[3][mu]) outward directed with central
    // pixel k_phi = 0. this ensures that a photon with zero impact
    // parameter will be in the center of the field of view.
    trial[0] = 1.;
    trial[1] = 1.;
    trial[2] = 0.;
    trial[3] = 0.;
    lower(trial, Gcon, Kcon);

    // set the y camera direction to be parallel to the projected
    // spin axis of the black hole (on the image plane defined to
    // be normal to the Kcon vector above).
    trial[0] = 0.;
    trial[1] = 0.;
    trial[2] = 1.;
    trial[3] = 0.;

    make_plasma_tetrad(Ucam, Kcon, trial, Gcov, Econ, Ecov);
  }

  // Construct outgoing wavevectors
  // xoff: allow arbitrary offset for e.g. ML training imgs
  // +0.5: project geodesics from px centers
  // -0.01: Prevent nasty artifact at 0 spin/phicam
  double rotcam = 0.;
  double dxoff = (di+i+0.5-0.01)/NX - 0.5;
  double dyoff = (dj+j+0.5)/NY - 0.5;
  Kcon_tetrad[0] = 0.;
  Kcon_tetrad[1] = (dxoff*cos(rotcam) - dyoff*sin(rotcam)) * fovx;
  Kcon_tetrad[2] = (dxoff*sin(rotcam) + dyoff*cos(rotcam)) * fovy;
  Kcon_tetrad[3] = 1.;

  if (Kcon_tetrad[0]) 
    ;

  // normalize kcon (in tetrad basis)
  double norm = sqrt(Kcon_tetrad[1]*Kcon_tetrad[1] + Kcon_tetrad[2]*Kcon_tetrad[2] + Kcon_tetrad[3]*Kcon_tetrad[3]);
  Kcon_tetrad[0] = 1.;
  Kcon_tetrad[1] /= norm;
  Kcon_tetrad[2] /= norm;
  Kcon_tetrad[3] /= norm;

  // translate into coordinate frame
  tetrad_to_coordinate(Econ, Kcon_tetrad, Kcon);

  // lower
  double Gcov[4][4];
  set_gcov(Xcam, Gcov);
  for (int i=0; i<4; i++) {
    Kcov[i] = 0.;
    for (int j=0; j<4; j++) {
      Kcov[i] += Gcov[i][j] * Kcon[j];
    }
  }

  // update camera position
  for (int mu=0; mu<4; ++mu) {
    X[mu] = Xcam[mu];
  }
}

void push_photon_hidden (double X[4], double Kcon[4])
{
  double Xhalf[4], Kconhalf[4];
  double dl = stepsize(X, Kcon);
  push_photon(X, Kcon, -dl, Xhalf, Kconhalf);
}

int push_photon_hidden_n_dl (double Xarr[][5], double Kcon[4], int n)
{  
  double X[4], Xhalf[4], Kconhalf[4];
  int k = 0;

  for (int i=0; i<4; ++i) X[i] = Xarr[0][i];

  while ( k < n && valid_geodesic(X,Kcon) ) {
    double dl = stepsize(X, Kcon);
    push_photon(X, Kcon, -dl, Xhalf, Kconhalf);
    for (int i=0; i<4; ++i) Xarr[k][i] = X[i];
    Xarr[k][4] = dl;
    k++;
  }

  return k;
}

int push_photon_hidden_n (double Xarr[][4], double Kcon[4], int n)
{
  double X[4], Xhalf[4], Kconhalf[4];
  int k = 0;

  for (int i=0; i<4; ++i) X[i] = Xarr[0][i];

  while ( k < n && valid_geodesic(X,Kcon) ) {
    double dl = stepsize(X, Kcon);
    push_photon(X, Kcon, -dl, Xhalf, Kconhalf);
    for (int i=0; i<4; ++i) Xarr[k][i] = X[i];
    k++;
  }

  return k;
}

void set_dl (double Xarr[][4], double dls[], int n) 
{
  // WARNING! This expects Xarr to be given in KS coordinates.
  double gcov[4][4];
  double dl;
  for (int i=1; i<n; ++i) {
    gcov_ks(Xarr[i][1], Xarr[i][2], gcov);
    dl = 0.;
    for (int mu=0; mu<4; ++mu) {
      dl += gcov[0][mu] * (Xarr[i-1][mu] - Xarr[i][mu]);
    }
    dls[i] = dl;
  }
}

void set_dl_native_coordinates (double Xarr[][4], double dls[], int n) 
{
  double gcov[4][4];
  double dl;
  for (int i=1; i<n; ++i) {
    set_gcov(Xarr[i], gcov);
    dl = 0.;
    for (int mu=0; mu<4; ++mu) {
      dl += gcov[0][mu] * (Xarr[i-1][mu] - Xarr[i][mu]);
    }
    dls[i] = dl;
  }
}

double solve_temp(double Ii, double ji, double jf, double dl)
{
  double javg = (ji + jf) / 2.;
  return Ii + javg * dl;
}
