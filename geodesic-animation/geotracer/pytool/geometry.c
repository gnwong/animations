
#include "geometry.h"
#include "mathfunc.h"
#include "model.h"

double eps = 0.01;

void set_eps(double e) {
  eps = e;
}

// TODO rewrite description
// generic connection routine, using numerical derivatives 
// Sets the spatial discretization in numerical derivatives: 
void get_connection(double X[4], double conn[4][4][4]) {

  double DEL = 1.e-7;

  double Xh[4], Xl[4];
  double gcon[4][4], gcov[4][4];
  double gh[4][4], gl[4][4];
  double tmp[4][4][4];

  set_gcov(X, gcov);
  invert_44(gcov, gcon);

  for (int k=0; k<4; ++k) {
    for (int l=0; l<4; ++l) {
      Xh[l] = X[l];
      Xl[l] = X[l];
    }
    Xh[k] += DEL;
    Xl[k] -= DEL;
    set_gcov(Xh, gh);
    set_gcov(Xl, gl);
    for (int i=0; i<4; ++i) {
      for (int j=0; j<4; ++j) {
        conn[i][j][k] =  (gh[i][j] - gl[i][j]) / (Xh[k] - Xl[k]);
      }
    }
  }

  // rearrange
  for (int i=0; i<4; ++i) {
    for (int j=0; j<4; ++j) {
      for (int k=0; k<4; ++k) {
        tmp[i][j][k] =  0.5 * (conn[j][i][k] + conn[k][i][j] - conn[k][j][i]);
      }
    }
  }

  // and raise index
  for (int i=0; i<4; ++i) {
    for (int j=0; j<4; ++j) {
      for (int k=0; k<4; ++k) {
        conn[i][j][k] = 0.;
        for (int l=0; l<4; ++l) {
          conn[i][j][k] += gcon[i][l] * tmp[l][j][k];
        }
      }
    }
  }
}

// 2nd order integrator. 
void push_photon(double X[4], double Kcon[4], double dl,double Xhalf[4],double Kconhalf[4]) {

  double lconn[4][4][4];
  double dKcon[4], Xh[4], Kconh[4];

  // take half step
  get_connection(X, lconn);
  for (int k=0; k<4; ++k) dKcon[k] = 0.;
  for (int i=0; i<4; ++i) {
    for (int j=0; j<4; ++j) {
      for (int k=0; k<4; ++k) {
        dKcon[k] -= 0.5 * dl * lconn[k][i][j] * Kcon[i] * Kcon[j];
      }
    }
  }

  // update and save
  for (int i=0; i<4; ++i) Kconh[i] = Kcon[i] + dKcon[i];
  for (int i=0; i<4; ++i) Xh[i] = X[i] + 0.5 * dl * Kcon[i];
  for (int i=0; i<4; ++i) Xhalf[i] = Xh[i];
  for (int i=0; i<4; ++i) Kconhalf[i] = Kconh[i];

  // finish out with full step
  get_connection(Xh, lconn);
  for (int k=0; k<4; ++k) dKcon[k] = 0.;
  for (int i=0; i<4; ++i) {
    for (int j=0; j<4; ++j) {
      for (int k=0; k<4; ++k) {
        dKcon[k] -= dl * lconn[k][i][j] * Kconh[i] * Kconh[j];

      }
    }
  }

  // update
  for (int i=0; i<4; ++i) Kcon[i] += dKcon[i];
  for (int i=0; i<4; ++i) X[i] += dl * Kconh[i];
}


#define MIN(A,B) (A<B?A:B)
double stepsize (double X[4], double Kcon[4]) {

  // these should be optimized away by the compiler, so we can write them
  // here to avoid muddying the namespace
  double small = 1.e-40; 

  // no need to touch anything below this line.
  double dlx1, dlx2, dlx3;
  double idlx1, idlx2, idlx3;

  dlx1 = eps / (fabs(Kcon[1]) + small*small);
  dlx2 = eps * MIN(X[2], 1. - X[2]) / (fabs(Kcon[2]) + small*small);
  dlx3 = eps / (fabs(Kcon[3]) + small*small);

  idlx1 = 1./(fabs(dlx1) + small*small);
  idlx2 = 1./(fabs(dlx2) + small*small);
  idlx3 = 1./(fabs(dlx3) + small*small);

  return 1. / (idlx1 + idlx2 + idlx3);
}
#undef MIN

/* input and vectors are contravariant (index up) */
void tetrad_to_coordinate(double Econ[4][4], double K_tetrad[4], double K[4])
{
    int l;

    for (l = 0; l < 4; l++) {
  K[l] = Econ[0][l] * K_tetrad[0] 
       + Econ[1][l] * K_tetrad[1] 
       + Econ[2][l] * K_tetrad[2] 
       + Econ[3][l] * K_tetrad[3];
    }

    return;
}


/* assumes gcov has been set first; returns sqrt{|g|} */
double gdet_func(double gcov[][4])
{

    int i,j;
    int permute[4];
    double gcovtmp[4][4];
    double gdet;
    int LU_decompose(double A[][4], int permute[]);

    for (i = 0; i < 4; i++) 
    for (j = 0; j < 4; j++) {
  gcovtmp[i][j] = gcov[i][j];
    }
    if (LU_decompose(gcovtmp, permute) != 0) {
  fprintf(stderr, "gdet_func(): singular matrix encountered! \n");
      //exit(1);
    }
    gdet = 1.;

    for (i = 0; i < 4; i++)
  gdet *= gcovtmp[i][i];

    return (sqrt(fabs(gdet)));
}


void lower(double ucon[4], double Gcov[4][4], double ucov[4]) {
    ucov[0] = Gcov[0][0] * ucon[0] + Gcov[0][1] * ucon[1]
            + Gcov[0][2] * ucon[2] + Gcov[0][3] * ucon[3];
    ucov[1] = Gcov[1][0] * ucon[0] + Gcov[1][1] * ucon[1]
            + Gcov[1][2] * ucon[2] + Gcov[1][3] * ucon[3];
    ucov[2] = Gcov[2][0] * ucon[0] + Gcov[2][1] * ucon[1]
            + Gcov[2][2] * ucon[2] + Gcov[2][3] * ucon[3];
    ucov[3] = Gcov[3][0] * ucon[0] + Gcov[3][1] * ucon[1]
            + Gcov[3][2] * ucon[2] + Gcov[3][3] * ucon[3];
}




void normalize(double vcon[4], double Gcov[4][4]) {
    int k, l;
    double norm;

    norm = 0.;
    for (k = 0; k < 4; k++)
  for (l = 0; l < 4; l++)
      norm += vcon[k] * vcon[l] * Gcov[k][l];

    norm = sqrt(fabs(norm));
    for (k = 0; k < 4; k++)
  vcon[k] /= norm;

    return;
}


#define SMALL_VECTOR (1.e-30)
void set_Econ_from_trial(double Econ[4], int defdir, double trial[4])
{
    double norm = 0.;
    int k;

    for (k = 0; k < 4; k++) 
        norm += fabs(trial[k]);
    for (k = 0; k < 4; k++) /* trial vector */
        if (norm <= SMALL_VECTOR) /* bad trial vector; default to defdir */
            Econ[k] = (k==defdir)?1:0;
        else
            Econ[k] = trial[k];

}


double check_handedness(double Econ[4][4], double Gcov[4][4])
{
    int i, j, k, l;

    double g = gdet_func(Gcov);

    /* check handedness */
    double dot = 0.;
    for (i = 0; i < 4; i++)
    for (j = 0; j < 4; j++)
        for (l = 0; l < 4; l++)
        for (k = 0; k < 4; k++) {
            dot += g * (j-i)*(k-i)*(l-i)*(k-j)*(l-j)*(l-k)/12 * // levi civita
            Econ[0][i] * Econ[1][j] * Econ[2][k] * Econ[3][l];
        }

    return (dot);
}

void project_out(double vcona[4], double vconb[4], double Gcov[4][4])
{

    double adotb, vconb_sq;
    int k, l;

    vconb_sq = 0.;
    for (k = 0; k < 4; k++)
  for (l = 0; l < 4; l++)
      vconb_sq += vconb[k] * vconb[l] * Gcov[k][l];

    adotb = 0.;
    for (k = 0; k < 4; k++)
  for (l = 0; l < 4; l++)
      adotb += vcona[k] * vconb[l] * Gcov[k][l];

    for (k = 0; k < 4; k++)
  vcona[k] -= vconb[k] * adotb / vconb_sq;

    return;
}


#define METRIC_eKS (1)
void make_plasma_tetrad(double Ucon[4], double Kcon[4],
      double Bcon[4], double Gcov[4][4],
      double Econ[4][4], double Ecov[4][4])
{
    int k, l;
    void normalize(double *vcon, double Gcov[4][4]);
    void project_out(double *vcona, double *vconb, double Gcov[4][4]);
    double check_handedness(double Econ[4][4],
          double Gcov[4][4]);

    void set_Econ_from_trial(double Econ[4], int defdir, double trial[4]);

    // start with time component along Ucon
    set_Econ_from_trial(Econ[0], 0, Ucon);
    normalize(Econ[0], Gcov);

    // use trial vector Kcon to set Econ[3] in usual way
    set_Econ_from_trial(Econ[3], 3, Kcon);
    project_out(Econ[3], Econ[0], Gcov);
    normalize(Econ[3], Gcov);

    // do usual thing for Econ[2] with Bcon
    set_Econ_from_trial(Econ[2], 2, Bcon);
    project_out(Econ[2], Econ[3], Gcov);
    normalize(Econ[2], Gcov);
    
    // set whatever's left to Econ[1]
    for (k = 0; k < 4; k++) Econ[1][k] = 1.;
    project_out(Econ[1], Econ[0], Gcov);
    project_out(Econ[1], Econ[2], Gcov);
    project_out(Econ[1], Econ[3], Gcov);
    normalize(Econ[1], Gcov);

    // check handedness
    double dot = check_handedness(Econ, Gcov);

    // less restrictive condition on geometry for eKS coordinates which are
    // used when the exotic is expected.
    if ((fabs(fabs(dot) - 1.) > 1.e-10 && METRIC_eKS == 0) ||
        (fabs(fabs(dot) - 1.) > 1.e-7  && METRIC_eKS == 1)) {
      fprintf(stderr, "that's odd: %g\n", fabs(dot) - 1.);
      fprintf(stderr, "Ucon[] = %e %e %e %e\n", Ucon[0], Ucon[1],Ucon[2],Ucon[3]);
      fprintf(stderr, "Kcon[] = %e %e %e %e\n", Kcon[0], Kcon[1],Kcon[2],Kcon[3]);
      fprintf(stderr, "Bcon[] = %e %e %e %e\n", Bcon[0], Bcon[1], Bcon[2], Bcon[3]);

      double ucov[4];
      lower(Ucon, Gcov, ucov);
      double udotu=0., udotb=0.;
      for (int mu=0; mu<4; ++mu) {
        udotu += Ucon[mu]*ucov[mu];
        udotb += Bcon[mu]*ucov[mu];
      }
      fprintf(stderr, "u.u = %g  u.b = %g\n", udotu, udotb);
    }

    // ensure handedness == right
    if (dot < 0.) {
      for (k = 0; k < 4; k++) Econ[1][k] *= -1.;
    }

    // make covariant
    for (k = 0; k < 4; k++) lower(Econ[k], Gcov, Ecov[k]);

    // raise tetrad basis index
    for (l = 0; l < 4; l++) Ecov[0][l] *= -1.;
}

