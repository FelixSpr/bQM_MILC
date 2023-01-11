// -----------------------------------------------------------------
// Update the momentum matrices
// Uncomment to print out debugging messages
//#define FORCE_DEBUG
#include "bQM_includes.h"
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Update mom with the bosonic force
double bosonic_force(Real eps) {
  register int i, j;
  register site *s;
  Real tr;
  double returnit = 0.0;
  matrix tmat, tmat2;
  msg_tag *tag[NSCALAR], *tag2[NSCALAR];
#ifdef DEBUG_CHECK
  anti_hermitmat tah;
#endif

  // Clear the force collectors
  FORALLSITES(i, s) {
    clear_mat(&(s->f_U));
    for (j = 0; j < NSCALAR; j++)
      clear_mat(&(s->f_X[j]));
  }

  // First we have the finite difference operator gauge derivative
  // Must transform as site variable so momenta can be exponentiated
  //   U(n) d/dU(n) Tr[2 U(t) X(t+1) Udag(t) X(t) - X(t+1) X(t+1) - X(t) X(t)]
  //     = 2 delta_{nt} U(n) X(t+1) Udag(t) X(t) = 2 U(n) X(n+1) Udag(n) X(n)
  for (j = 0; j < NSCALAR; j++) {
    tag[j] = start_gather_site(F_OFFSET(X[j]), sizeof(matrix),
                               TUP, EVENANDODD, gen_pt[j]);
  }

  for (j = 0; j < NSCALAR; j++) {
    // For scalar force term, compute and gather Udag(n-1) X(n-1) U(n-1)
    FORALLSITES(i, s) {
      mult_nn(&(s->X[j]), &(s->link), &tmat);
      mult_an(&(s->link), &tmat, &(temp_X[j][i]));
    }
    tag2[j] = start_gather_field(temp_X[j], sizeof(matrix),
                                 TDOWN, EVENANDODD, gen_pt[NSCALAR + j]);
  }

  for (j = 0; j < NSCALAR; j++) {   // X(n+1) = gen_pt[j]
    wait_gather(tag[j]);
    FORALLSITES(i, s) {
      mult_na((matrix *)(gen_pt[j][i]), &(s->link), &tmat);
      mult_nn(&(s->link), &tmat, &tmat2);
      mult_nn_sum(&(s->X[j]), &tmat2, &(s->f_U));
    }
  }

  // Take adjoint and update the gauge momenta
  // Make them anti-hermitian
  // Include overall factor of 2
  // !!! Another factor of 2 needed for conservation (real vs. complex?)
  // Compute average gauge force in same loop
  tr = 4.0 * eps;
  FORALLSITES(i, s) {
    uncompress_anti_hermitian(&(s->mom), &tmat);
    scalar_mult_dif_matrix(&(s->f_U), tr, &tmat);
    make_anti_hermitian(&tmat, &(s->mom));
    returnit += realtrace(&(s->f_U), &(s->f_U));
  }

  // This is the finite difference operator scalar derivative
  //   d/dX(n) Tr[X(t) U(t) X(t+1) Udag(t) + X(t+1) Udag(t) X(t) U(t)
  //              - X(t+1) X(t+1) - X(t) X(t)]
  //     = 2 delta_{nt} U(t) X(t+1) Udag(t)
  //       + 2 delta_{n(t+1)} Udag(t) X(t) U(t)
  //       - 2 delta_{n(t+1)} X(t+1) - 2 delta_{nt} X(t)
  //     = 2 [U(n) X(n+1) Udag(n) + Udag(n-1) X(n-1) U(n-1) - 2 X(n)]
  for (j = 0; j < NSCALAR; j++) {
    wait_gather(tag2[j]);
    FORALLSITES(i, s) {
      // Initialize force with on-site -2X(n)
      scalar_mult_matrix(&(s->X[j]), -2.0, &(s->f_X[j]));

      // Add forward hopping term using X(n+1) = gen_pt[j]
      mult_na((matrix *)(gen_pt[j][i]), &(s->link), &tmat);
      mult_nn_sum(&(s->link), &tmat, &(s->f_X[j]));

      // Add backward hopping term
      //   Udag(n-1) X(n-1) U(n-1) = gen_pt[NSCALAR + j]
      sum_matrix((matrix *)(gen_pt[NSCALAR + j][i]), &(s->f_X[j]));

      scalar_mult_matrix(&(s->f_X[j]), 2.0, &(s->f_X[j]));
    }
    cleanup_gather(tag[j]);
    cleanup_gather(tag2[j]);
  }

  // The pure scalar stuff: simple d/dX_i(n) -X_j(t)^2 = -2 X_i(n)
  // Absorb factor of two into tr = 2 * (2 + omega^2)
  tr = 2.0 * (2.0 + omega * omega);
  FORALLSITES(i, s) {
    for (j = 0; j < NSCALAR; j++)
      scalar_mult_dif_matrix(&(s->X[j]), tr, &(s->f_X[j]));
  }

  // Take adjoint and update the scalar momenta
  // Subtract to reproduce -Adj(f_X)
  // Compute average scalar force in same loop (combine with gauge from above)
  FORALLSITES(i, s) {
    for (j = 0; j < NSCALAR; j++) {
#ifdef DEBUG_CHECK
      // Make f_X traceless anti-hermitian, which it should be already
      make_anti_hermitian(&(s->f_X[j]), &tah);
      uncompress_anti_hermitian(&tah, &(s->f_X[j]));
#endif
      scalar_mult_sum_matrix(&(s->f_X[j]), eps, &(s->mom_X[j]));
      returnit += realtrace(&(s->f_X[j]), &(s->f_X[j]));
    }
  }
  g_doublesum(&returnit);

  return (eps * sqrt(returnit) / (double)nt);
}
// -----------------------------------------------------------------
