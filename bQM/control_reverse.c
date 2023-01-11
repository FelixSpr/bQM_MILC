// -----------------------------------------------------------------
// Main procedure for bosonic quantum mechanics reversibility test
#define CONTROL
#include "bQM_includes.h"

// A bit hacky... this is defined in update_o.c, which we will link
int update_step();
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// A bit hacky... copy update() without discarding pseudofermions at the end
// Can also strip out some HMC-related stuff since we always accept
void persistent_update() {
  int n;
  double startaction, endaction, change;

  // Refresh the momenta
  ranmom();

  // Find initial action and do microcanonical updating
  // No need to save links, since will automatically accept
  startaction = action(source, psim);
  bnorm = 0.0;
  max_bf = 0.0;
  update_step(source, psim);

  // Find new action, always accept
  endaction = action(source, psim);
  change = endaction - startaction;

  // Warn about overflow
  if (fabs((double)change) > 1e20) {
    node0_printf("WARNING: Correcting apparent overflow: Delta S = %.4g\n",
                 change);
    change = 1e20;
  }
  node0_printf("delta S = %.4g, start S = %.12g, end S = %.12g\n",
               change, startaction, endaction);

  if (traj_length > 0) {
    node0_printf("MONITOR_FORCE %.4g %.4g\n",
                 bnorm / (double)(2.0 * nsteps), max_bf);
  }
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// A copy of persistent_update(), reversing instead of refreshing the momenta
// and not generating a new pseudofermion configuration
void reverse() {
  register int i, j, n;
  register site *s;
  double startaction, endaction, change;

  // Reverse the momenta (anti_hermitmat defined in include/bQM.h)
  FORALLSITES(i, s) {
    for (j = 0; j < NCOL; j++)
      s->mom.im_diag[j] *= -1.0;
    for (j = 0; j < N_OFFDIAG; j++)
      CNEGATE(s->mom.m[j], s->mom.m[j]);

    for (j = 0; j < NSCALAR; j++)
      scalar_mult_matrix(&(s->mom_X[j]), -1.0, &(s->mom_X[j]));
  }

  // Find initial action and do microcanonical updating
  // No need to save links, since will automatically accept
  startaction = action(source, psim);
  bnorm = 0.0;
  max_bf = 0.0;
  update_step(source, psim);

  // Find new action, always accept
  endaction = action(source, psim);
  change = endaction - startaction;

  // Warn about overflow
  if (fabs((double)change) > 1e20) {
    node0_printf("WARNING: Correcting apparent overflow: Delta S = %.4g\n",
                 change);
    change = 1e20;
  }
  node0_printf("delta S = %.4g, start S = %.12g, end S = %.12g\n",
               change, startaction, endaction);

  if (traj_length > 0) {
    node0_printf("MONITOR_FORCE %.4g %.4g\n",
                 bnorm / (double)(2.0 * nsteps), max_bf);
  }
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
int main(int argc, char *argv[]) {
  int prompt, j;
  Real eps;
  double start_act, end_act, dtime, Xtr[NSCALAR], Xtr_ave, Xtr_width;
  double ave_eigs[NCOL], eig_widths[NCOL], min_eigs[NCOL], max_eigs[NCOL];
  complex plp = cmplx(99.0, 99.0);

  // Setup
  setlinebuf(stdout); // DEBUG
  initialize_machine(&argc, &argv);
  // Remap standard I/O
  if (remap_stdio_from_args(argc, argv) == 1)
    terminate(1);

  g_sync();
  prompt = setup();
  setup_lambda();
  setup_gamma();
  setup_rhmc();

  // Always accept
#ifdef HMC_ALGORITHM
#undef HMC_ALGORITHM
#endif

  // Load input and run
  if (readin(prompt) != 0) {
    node0_printf("ERROR in readin, aborting\n");
    terminate(1);
  }
  dtime = -dclock();

  // Bosonic observables at start of trajectory
  // Action, scalar squares, scalar eigenvalues and Polyakov loop
  start_act = bosonic_action(&(Xtr[0]), &(Xtr[1]), &(Xtr[2]), &(Xtr[3]));
  node0_printf("ACT %.8g %.8g %.8g %.8g\n",
               start_act / (double)nt, Xtr[0] / (double)nt,
               Xtr[1] / (double)nt, Xtr[2] / (double)nt);

  Xtr_ave = scalar_trace(Xtr, &Xtr_width);
  node0_printf("SCALAR SQUARES");
  for (j = 0; j < NSCALAR; j++)
    node0_printf(" %.6g", Xtr[j]);
  node0_printf(" %.6g %.6g\n", Xtr_ave, Xtr_width);

  scalar_eig(ave_eigs, eig_widths, min_eigs, max_eigs);
  for (j = 0; j < NCOL; j++) {
    node0_printf("SCALAR_EIG %d %.6g %.6g %.6g %.6g\n",
                 j, ave_eigs[j], eig_widths[j], min_eigs[j], max_eigs[j]);
  }

  plp = ploop_eig();
  node0_printf("POLYA %.8g %.8g\n", plp.real, plp.imag);

  // Evolve forward for one trajectory
  eps = traj_length / (Real)nsteps;
  node0_printf("eps %.4g\n", eps);
  persistent_update(source, psim);

  // Bosonic observables at end of trajectory
  end_act = bosonic_action(&(Xtr[0]), &(Xtr[1]), &(Xtr[2]), &(Xtr[3]));
  node0_printf("ACT %.8g %.8g %.8g %.8g\n",
               end_act / (double)nt, Xtr[0] / (double)nt,
               Xtr[1] / (double)nt, Xtr[2] / (double)nt);

  Xtr_ave = scalar_trace(Xtr, &Xtr_width);
  node0_printf("SCALAR SQUARES");
  for (j = 0; j < NSCALAR; j++)
    node0_printf(" %.6g", Xtr[j]);
  node0_printf(" %.6g %.6g\n", Xtr_ave, Xtr_width);

  scalar_eig(ave_eigs, eig_widths, min_eigs, max_eigs);
  for (j = 0; j < NCOL; j++) {
    node0_printf("SCALAR_EIG %d %.6g %.6g %.6g %.6g\n",
                 j, ave_eigs[j], eig_widths[j], min_eigs[j], max_eigs[j]);
  }

  plp = ploop_eig();
  node0_printf("POLYA %.8g %.8g\n", plp.real, plp.imag);

  // Reverse momenta and evolve backwards for one trajectory
  reverse(source, psim);

  // Bosonic observables hopefully back at the start of the trajectory
  end_act = bosonic_action(&(Xtr[0]), &(Xtr[1]), &(Xtr[2]), &(Xtr[3]));
  node0_printf("ACT %.8g %.8g %.8g %.8g\n",
               end_act / (double)nt, Xtr[0] / (double)nt,
               Xtr[1] / (double)nt, Xtr[2] / (double)nt);

  Xtr_ave = scalar_trace(Xtr, &Xtr_width);
  node0_printf("SCALAR SQUARES");
  for (j = 0; j < NSCALAR; j++)
    node0_printf(" %.6g", Xtr[j]);
  node0_printf(" %.6g %.6g\n", Xtr_ave, Xtr_width);

  scalar_eig(ave_eigs, eig_widths, min_eigs, max_eigs);
  for (j = 0; j < NCOL; j++) {
    node0_printf("SCALAR_EIG %d %.6g %.6g %.6g %.6g\n",
                 j, ave_eigs[j], eig_widths[j], min_eigs[j], max_eigs[j]);
  }

  plp = ploop_eig();
  node0_printf("POLYA %.8g %.8g\n", plp.real, plp.imag);

  // Done
  node0_printf("RUNNING COMPLETED\n");
  node0_printf("Initial action: %.8g\n", start_act);
  node0_printf("Final action:   %.8g\n", end_act);
  node0_printf("Difference:     %.4g\n", fabs(end_act - start_act));
  dtime += dclock();
  node0_printf("Time = %.4g seconds\n", dtime);
  fflush(stdout);
  return 0;
}
// -----------------------------------------------------------------
