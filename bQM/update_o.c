// -----------------------------------------------------------------
// Update lattice
// Omelyan integrator multiscale following CPC 174:87 (2006)

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
// Omelyan version; ``dirty'' speeded-up version
double update_bosonic_step(Real eps) {
  int n = nsteps, i;
  double norm;
#ifdef UPDATE_DEBUG
  double td, td2;
#endif

#ifdef UPDATE_DEBUG
  node0_printf("gauge %d steps %.4g dt\n", n, eps);
#endif
  norm = bosonic_force(eps * LAMBDA);
  for (i = 1; i <= n; i++) {
    update_u(0.5 * eps);
    norm += bosonic_force(eps * LAMBDA_MID);
    update_u(0.5 * eps);
    if (i < n)
      norm += bosonic_force(eps * TWO_LAMBDA);

    else
      norm += bosonic_force(eps * LAMBDA);
  }

  // Reunitarize the gauge field and re-anti-hermitianize the scalars
#ifdef UPDATE_DEBUG
  td = check_unitarity();
  g_floatmax(&td);
#endif
  reunitarize();
  reantihermize();
#ifdef UPDATE_DEBUG
  td2 = check_unitarity();
  g_floatmax(&td2);
  node0_printf("Reunitarized after boson update step.  ");
  node0_printf("Max deviation %.2g changed to %.2g\n", td, td2);
#endif

  return (norm / n);
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
void update_step() {
  int i_multi0;
  Real eps, tr;

  eps = traj_length / (Real)nsteps;

  for (i_multi0 = 1; i_multi0 <= nsteps; i_multi0++) {
    tr = update_bosonic_step(eps);
#ifdef UPDATE_DEBUG
    node0_printf("Step %d - 1 Action %.4g Force %.4g\n", i_multi0,
                 action(), tr);
#endif
    bnorm += tr;
    if (tr > max_bf)
      max_bf = tr;

    tr = update_bosonic_step(eps);
#ifdef UPDATE_DEBUG
    node0_printf("Step %d - 2 Action %.4g Force %.4g\n", i_multi0,
                 action(), tr);
#endif
    bnorm += tr;
    if (tr > max_bf)
      max_bf = tr;

  }
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
void update() {
  int n;
  double startaction, endaction, change;

  // Refresh the momenta
  ranmom();

  // Find initial action
  startaction = action();
  bnorm = 0.0;
  max_bf = 0.0;

#ifdef HMC_ALGORITHM
  Real xrandom;   // For accept/reject test
  // Copy link field and scalars to old_link and old_X
  copy_bosons(PLUS);
#endif
  // Do microcanonical updating
  update_step();

  // Find ending action
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
    if (traj_length > 0.0) {
      // Restore link field and scalars from old_link and old_X
      copy_bosons(MINUS);
    }
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
                 bnorm / (double)(2.0 * nsteps), max_bf);
  }
}
// -----------------------------------------------------------------
