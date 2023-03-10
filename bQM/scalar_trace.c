// -----------------------------------------------------------------
// Measure the average value of Tr[X[i] X[i]] / N
// as well as the width sqrt(<tr^2> - <tr>^2) of its distribution
#include "bQM_includes.h"

double scalar_trace(double *Xtr, double *Xwidth) {
  register int i, j;
  register site *s;
  double Xtr_ave = 0.0, XtrSq = 0.0, td;

  for (j = 0; j < NSCALAR; j++) {
    Xtr[j] = 0.0;
    FORALLSITES(i, s) {
      // Take adjoint of first to get rid of overall negative sign
      td = realtrace(&(s->X[j]), &(s->X[j]));
      Xtr[j] += td;
      XtrSq += td * td;
    }
    Xtr[j] *= one_ov_N / ((double)nt);
    g_doublesum(&(Xtr[j]));
    Xtr_ave += Xtr[j];
  }
  Xtr_ave /= (double)NSCALAR;
  XtrSq *= one_ov_N * one_ov_N / ((double)nt * NSCALAR);
  g_doublesum(&XtrSq);
  *Xwidth = sqrt(XtrSq - Xtr_ave * Xtr_ave);

  return Xtr_ave;
}
// -----------------------------------------------------------------
