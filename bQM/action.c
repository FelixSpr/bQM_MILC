// -----------------------------------------------------------------
// Measure total action, as needed by the hybrid Monte Carlo algorithm
#include "bQM_includes.h"
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Bosonic contribution to the action
double bosonic_action() {
  register int i,l;
  register site *s;
  int j, k;
  double b_action = 0.0;
  double sqterms = 0.0;
  matrix tmat, tmat2;
  msg_tag *tag[NSCALAR];

  // Scalar kinetic term -Tr[D_t X(t)]^2
  //   -Tr[U(t) X(t+1) Udag(t) - X(t)]^2
  //     = Tr[2 X(t) U(t) X(t+1) Udag(t) - X(t+1) X(t+1) - X(t) X(t)]
  // Sum over t --> 2 Tr[Udag(t) X(t) U(t) X(t+1) - X(t) X(t)]
  for (j = 0; j < NSCALAR; j++) {
    tag[j] = start_gather_site(F_OFFSET(X[j]), sizeof(matrix),
                               TUP, EVENANDODD, gen_pt[j]);
  }

  // On-site piece of scalar kinetic term
  // (Has same form as some scalar potential terms, so will re-use below)
  FORALLSITES(i, s) {
    for (j = 0; j < NSCALAR; j++)
      sqterms -= (double)realtrace_nn(&(s->X[j]), &(s->X[j]));
  }
  b_action = (2.0+omega*omega)*sqterms;

  // Nearest-neighbor piece of scalar kinetic term
  for (j = 0; j < NSCALAR; j++) {
    wait_gather(tag[j]);
    FORALLSITES(i, s) {

      mult_nn(&(s->link), (matrix *)(gen_pt[j][i]), &tmat);
      mult_na(&tmat, &(s->link), &tmat2);
      b_action += 2.0*(double)realtrace_nn(&(s->X[j]), &tmat2);
    }
    cleanup_gather(tag[j]);
  }
  

  // Scalar potential terms
  // TODO: REPLACE WITH BQM POTENTIAL..
  //b_action += ...;

  g_doublesum(&b_action);
  return b_action;
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Gauge and scalar momenta contribution to the action
// Helper routine computes agnitude squared of an anti-hermition matrix
// including the factor of 1/2 in the effective hamiltonian
Real ahmat_mag_sq(anti_hermitmat *ah) {
  register int i;
  register Real sum;

  sum = ah->im_diag[0] * ah->im_diag[0];
  for (i = 1; i < NCOL; i++)
    sum += ah->im_diag[i] * ah->im_diag[i];
  sum *= 0.5;

  for (i = 0; i < N_OFFDIAG; i++) {
    sum += ah->m[i].real * ah->m[i].real;
    sum += ah->m[i].imag * ah->m[i].imag;
  }

  return sum;
}

double gauge_mom_action() {
  register int i;
  register site *s;
  double sum = 0.0;

  FORALLSITES(i, s)
    sum += (double)ahmat_mag_sq(&(s->mom));

  g_doublesum(&sum);
  return sum;
}

double scalar_mom_action() {
  register int i, j;
  register site *s;
  double sum = 0.0;

  FORALLSITES(i, s) {
    for (j = 0; j < NSCALAR; j++)
      sum += (double)realtrace(&(s->mom_X[j]), &(s->mom_X[j]));
  }
  g_doublesum(&sum);
  return 0.5 * sum;
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Print out zeros for pieces of the action that aren't included
double action(matrix ***src, matrix ****sol) {
  double p_act, so3_act, so6_act, comm_act, Myers_act, total;

  // Includes so3, so6, Myers and kinetic
  total = bosonic_action(&so3_act, &so6_act, &comm_act, &Myers_act);
  node0_printf("action: so3 %.8g so6 %.8g comm %.8g Myers %.8g boson %.8g ",
               so3_act, so6_act, comm_act, Myers_act, total);

  p_act = gauge_mom_action();
  node0_printf("Umom %.8g ", p_act);
  total += p_act;
  p_act = scalar_mom_action();
  node0_printf("Xmom %.8g ", p_act);
  total += p_act;
  node0_printf("sum %.8g\n", total);
  return total;
}
// -----------------------------------------------------------------
