// -----------------------------------------------------------------
// Update lattice
// Leapfrog integrator
// Begin at "integral" time, with H and U evaluated at the same time

// Uncomment to print out debugging messages
//#define UPDATE_DEBUG
#include "bQM_includes.h"
#ifdef HAVE_IEEEFP_H
#include <ieeefp.h>         // For "finite"
#endif
// -----------------------------------------------------------------



// -----------------------------------------------------------------
void update_u(Real eps) {
  register int i, j;
  register site *s;
  register Real t2, t3, t4, t5, t6, t7, t8;
  matrix tmat, tmat2, tmp_mom;

  // Calculate newU = exp(p).U
  // Go to eighth order in the exponential:
  //   exp(p) * U = (1 + p + p^2/2 + p^3/6 ...) * U
  //              = U + p*(U + (p/2)*(U + (p/3)*( ... )))
  // Take divisions out of site loop (can't be done by compiler)
  t2 = eps / 2.0;
  t3 = eps / 3.0;
  t4 = eps / 4.0;
  t5 = eps / 5.0;
  t6 = eps / 6.0;
  t7 = eps / 7.0;
  t8 = eps / 8.0;

  FORALLSITES(i, s) {
    uncompress_anti_hermitian(&(s->mom), &tmp_mom);
    mult_nn(&tmp_mom, &(s->link), &tmat);
    scalar_mult_add_matrix(&(s->link), &tmat, t8, &tmat2);

    mult_nn(&tmp_mom, &tmat2, &tmat);
    scalar_mult_add_matrix(&(s->link), &tmat, t7, &tmat2);

    mult_nn(&tmp_mom, &tmat2, &tmat);
    scalar_mult_add_matrix(&(s->link), &tmat, t6, &tmat2);

    mult_nn(&tmp_mom, &tmat2, &tmat);
    scalar_mult_add_matrix(&(s->link), &tmat, t5, &tmat2);

    mult_nn(&tmp_mom, &tmat2, &tmat);
    scalar_mult_add_matrix(&(s->link), &tmat, t4, &tmat2);

    mult_nn(&tmp_mom, &tmat2, &tmat);
    scalar_mult_add_matrix(&(s->link), &tmat, t3, &tmat2);

    mult_nn(&tmp_mom, &tmat2, &tmat);
    scalar_mult_add_matrix(&(s->link), &tmat, t2, &tmat2);

    mult_nn(&tmp_mom, &tmat2, &tmat);
    scalar_mult_sum_matrix(&tmat, eps, &(s->link));

    for (j = 0; j < NSCALAR; j++)
      scalar_mult_sum_matrix(&(s->mom_X[j]), eps, &(s->X[j]));
  }
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
void update_step() {
  int step;
  Real eps = traj_length / (Real)nsteps, tr;
  node0_printf("eps %.4g\n", eps);

  // First u(t/2)
  update_u(0.5 * eps);
  
  for (step = 0; step < nsteps; step++) {
    //action();
    // Inner steps p(t) u(t)
    tr = bosonic_force(eps);
    bnorm += tr;
    if (tr > max_bf)
      max_bf = tr;

    if (step < nsteps - 1)
      update_u(eps);
    else                // Final u(t/2)
      update_u(0.5 * eps);
  }

  // Reunitarize the gauge field and re-anti-hermitianize the scalars
  reunitarize();
  reantihermize();
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
void update() {
  int i, j, k, n;
  site *s;
  double startaction, endaction, change;

  // Refresh the momenta
  ranmom();

  // Find initial action
  startaction = action();
  bnorm = 0.0;
  max_bf = 0.0;

#ifdef HMC_ALGORITHM
  Real xrandom;   // For accept/reject test
  // Copy link field to old_link
  copy_bosons(PLUS);
#endif
  // Do microcanonical updating
  update_step();

  // Find ending action
  // Since update_step ended on a gauge update,
  // need to do conjugate gradient to get (Mdag M)^(-1 / 4) chi
  endaction = action();
  change = endaction - startaction;
#ifdef HMC_ALGORITHM
  // Reject configurations giving overflow
#ifndef HAVE_IEEEFP_H
  if (fabs((double)change) > 1e20) {
#else
  if (!finite((double)change)) {
#endif
    node0_printf("WARNING: Correcting Apparent Overflow: Delta S = %.4g\n",
                 change);
    change = 1.0e20;
  }

  // Decide whether to accept, if not, copy old link field back
  // Careful -- must generate only one random number for whole lattice
  if (this_node == 0)
    xrandom = myrand(&node_prn);
  broadcast_float(&xrandom);
  if (exp(-change) < (double)xrandom) {
    if (traj_length > 0.0)
      copy_bosons(MINUS);
    node0_printf("REJECT: delta S = %.4g start S = %.12g end S = %.12g\n",
                 change, startaction, endaction);
  }
  else {
    node0_printf("ACCEPT: delta S = %.4g start S = %.12g end S = %.12g\n",
                 change, startaction, endaction);
  }
#else
  // Only print check if not doing HMC
  node0_printf("CHECK: delta S = %.4g\n", (double)(change));
#endif // ifdef HMC

  if (traj_length > 0) {
    node0_printf("MONITOR_FORCE %.4g %.4g\n",
                 bnorm / (double)(2 * nsteps), max_bf);
  }
}
// -----------------------------------------------------------------
