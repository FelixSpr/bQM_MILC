// -----------------------------------------------------------------
// Main procedure for bosonic quantum mechanics evolution and measurements
#define CONTROL
#include "bQM_includes.h"

int main(int argc, char *argv[]) {
  int prompt, j;
  int traj_done, Nmeas = 0;
  Real eps;
  double b_act, dtime, Xtr[NSCALAR], Xtr_ave, Xtr_width;
  double poloop_real=0;
  double poloop_imag=0;
  double poloop_abs=0;
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

  // Load input and run
  if (readin(prompt) != 0) {
    node0_printf("ERROR in readin, aborting\n");
    terminate(1);
  }
  dtime = -dclock();

  // Check: compute initial bosonic action and scalar squares
  //b_act = bosonic_action(&(Xtr[0]), &(Xtr[1]), &(Xtr[2]), &(Xtr[3]));
  b_act = bosonic_action();
  node0_printf("START %.8g\n", b_act / (double)nt);

  Xtr_ave = scalar_trace(Xtr, &Xtr_width);
  node0_printf("SCALAR SQUARES");
  for (j = 0; j < NSCALAR; j++)
    node0_printf(" %.6g", Xtr[j]);
  node0_printf(" %.6g %.6g\n", Xtr_ave, Xtr_width);

  // Perform warmup trajectories
  eps = traj_length / (Real)nsteps;
  node0_printf("eps %.4g\n", eps);
  for (traj_done = 0; traj_done < warms; traj_done++)
    update();
  node0_printf("WARMUPS COMPLETED\n");

  // Perform trajectories, reunitarizations and measurements
  for (traj_done = 0; traj_done < trajecs; traj_done++) {
    update();

    // Do "local" measurements every trajectory!
    // Tr[X^2] / N
    Xtr_ave = scalar_trace(Xtr, &Xtr_width);
    //node0_printf("SCALAR SQUARES");
    //for (j = 0; j < NSCALAR; j++)
    //  node0_printf(" %.6g", Xtr[j]);
    //node0_printf(" %.6g %.6g\n", Xtr_ave, Xtr_width);

    // Polyakov loop eigenvalues and trace
    // Format: GMES Re(Polyakov) Im(Poyakov)
    plp = ploop_eig();
    node0_printf("GMES %.8g %.8g\n", plp.real, plp.imag);
    poloop_real += plp.real;
    poloop_imag += plp.imag;
    poloop_abs += sqrt(pow(plp.real,2)+pow(plp.imag,2));
    node0_printf("Poloop RMS %.8g\n", sqrt(pow(plp.real,2)+pow(plp.imag,2)));
    // Bosonic action
    //b_act = bosonic_action(&(Xtr[0]), &(Xtr[1]), &(Xtr[2]), &(Xtr[3]));
    b_act = bosonic_action();
    node0_printf("b_act/nt %.8g\n",b_act / (double)nt);
    //node0_printf("b_act/nt Xtr[0]/nt Xtr[1]/nt Xtr[2]/nt %.8g %.8g %.8g %.8g\n",
    //             b_act / (double)nt, Xtr[0] / (double)nt,
    //             Xtr[1] / (double)nt, Xtr[2] / (double)nt);

    // Monitor scalar eigenvalues
    // Format: SCALAR_EIG # ave width min max
    //scalar_eig(ave_eigs, eig_widths, min_eigs, max_eigs);
    //for (j = 0; j < NCOL; j++) {
    //  node0_printf("SCALAR_EIG %d %.6g %.6g %.6g %.6g\n",
    //               j, ave_eigs[j], eig_widths[j], min_eigs[j], max_eigs[j]);
    //}

    // Less frequent measurements every "propinterval" trajectories
    // None at the moment
//    if ((traj_done % propinterval) == (propinterval - 1)) {
//    }
    fflush(stdout);
  }
  node0_printf("RUNNING COMPLETED\n");
  node0_printf("GMES %.8g %.8g\n", poloop_real/(double)trajecs, poloop_imag/(double)trajecs);
  node0_printf("Poloop RMS %.8g\n", poloop_abs/(double)trajecs);
  // Check: compute final bosonic action
  //b_act = bosonic_action(&(Xtr[0]), &(Xtr[1]), &(Xtr[2]), &(Xtr[3]));
  b_act = bosonic_action();
  node0_printf("STOP %.8g\n", b_act / (double)nt);
  dtime += dclock();
  node0_printf("\nTime = %.4g seconds\n", dtime);
  fflush(stdout);

  // Save lattice if requested
  if (saveflag != FORGET)
    save_lattice(saveflag, savefile);
  normal_exit(0);         // Needed by at least some clusters
  return 0;
}
// -----------------------------------------------------------------
