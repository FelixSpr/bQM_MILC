// -----------------------------------------------------------------
// Bosonic quantum mechanics setup
#include "bQM_includes.h"

#define IF_OK if (status == 0)

// Each node has a params structure for passing simulation parameters
#include "params.h"
params par_buf;
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// On node zero, read lattice size and seed, and send to others
int initial_set() {
  int prompt = 0, status = 0;
  if (mynode() == 0) {
    // Print banner
    printf("Bosonic QM, Nc = %d\n", NCOL);
    printf("Microcanonical simulation with refreshing\n");
    printf("Machine = %s, with %d nodes\n", machine_type(), numnodes());
#ifdef HMC_ALGORITHM
    printf("Hybrid Monte Carlo algorithm\n");
#endif
#ifdef PHI_ALGORITHM
    printf("Phi algorithm\n");
#else   // Quit!
    printf("Only works for phi algorithm\n");
    exit(1);
#endif
    time_stamp("start");
    status = get_prompt(stdin, &prompt);

    IF_OK status += get_i(stdin, prompt, "nt", &par_buf.nt);
    IF_OK status += get_i(stdin, prompt, "iseed", &par_buf.iseed);

    if (status > 0)
      par_buf.stopflag = 1;
    else
      par_buf.stopflag = 0;
  }

  // Broadcast parameter buffer from node 0 to all other nodes
  broadcast_bytes((char *)&par_buf, sizeof(par_buf));
  if (par_buf.stopflag != 0)
    normal_exit(0);

  nt = par_buf.nt;
  iseed = par_buf.iseed;

  // Lattice volume sanity checks, including dimensional reduction
  if (mynode() == 0) {
    if (nt < 1) {
      printf("nt must be positive\n");
      exit(1);
    }
  }

  this_node = mynode();
  number_of_nodes = numnodes();
  one_ov_N = 1.0 / (Real)NCOL;
  return prompt;
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Allocate space for fields
void make_fields() {
  // Temporary matrices and Fermions
  Real size = (Real)((2.0 + NSCALAR) * sizeof(matrix));
  FIELD_ALLOC(tempmat, matrix);
  FIELD_ALLOC(tempmat2, matrix);
  FIELD_ALLOC_VEC(temp_X, matrix, NSCALAR);

  size *= sites_on_node;
  node0_printf("Mallocing %.1f MBytes per core for fields\n", size / 1e6);
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
int setup() {
  int prompt;

  // Print banner, get volume and seed
  prompt = initial_set();
  // Initialize the node random number generator
  initialize_prn(&node_prn, iseed, nt + mynode());
  // Initialize the layout functions, which decide where sites live
  setup_layout();
  // Allocate space for lattice, set up coordinate fields
  make_lattice();
  // Set up neighbor pointers and comlink structures
  make_nn_gathers();

  // Allocate space for fields
  make_fields();

  return prompt;
}
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Read in parameters for Monte Carlo
// prompt=1 indicates prompts are to be given for input
int readin(int prompt) {
  int status;

  // On node zero, read parameters and send to all other nodes
  if (this_node == 0) {
    printf("\n\n");
    status = 0;

    // Warms, trajecs
    IF_OK status += get_i(stdin, prompt, "warms", &par_buf.warms);
    IF_OK status += get_i(stdin, prompt, "trajecs", &par_buf.trajecs);
    IF_OK status += get_f(stdin, prompt, "traj_length", &par_buf.traj_length);

    // Number of steps
    IF_OK status += get_i(stdin, prompt, "nstep", &par_buf.nsteps);

    // Trajectories between expensive measurements
    IF_OK status += get_i(stdin, prompt, "traj_between_meas",
                          &par_buf.propinterval);

    // beta, omega
    IF_OK status += get_f(stdin, prompt, "beta", &par_buf.beta);
    IF_OK status += get_f(stdin, prompt, "omega", &par_buf.omega);

    // Find out what kind of starting lattice to use
    IF_OK status += ask_starting_lattice(stdin, prompt, &par_buf.startflag,
                                         par_buf.startfile);

    // Find out what to do with lattice at end
    IF_OK status += ask_ending_lattice(stdin, prompt, &(par_buf.saveflag),
                                       par_buf.savefile);

    if (status > 0)
      par_buf.stopflag = 1;
    else
      par_buf.stopflag = 0;
  }

  // Broadcast parameter buffer from node0 to all other nodes
  broadcast_bytes((char *)&par_buf, sizeof(par_buf));
  if (par_buf.stopflag != 0)
    normal_exit(0);

  warms = par_buf.warms;
  trajecs = par_buf.trajecs;
  traj_length = par_buf.traj_length;
  nsteps = par_buf.nsteps;
  propinterval = par_buf.propinterval;

  beta = par_buf.beta;
  omega = par_buf.omega;

  startflag = par_buf.startflag;
  saveflag = par_buf.saveflag;
  strcpy(startfile, par_buf.startfile);
  strcpy(savefile, par_buf.savefile);

  // Allocate some more arrays to be used by LAPACK in scalar_eig.c
  // and in generic/reunitarize.c
  store = malloc(sizeof *store * 2 * NCOL * NCOL);
  eigs = malloc(sizeof *eigs * NCOL);
  work = malloc(sizeof *work * 4 * NCOL);
  Rwork = malloc(sizeof *Rwork * (3 * NCOL - 2));
  // Following are used only for reunitarization
  left = malloc(sizeof *left * 2 * NCOL * NCOL);
  right = malloc(sizeof *right * 2 * NCOL * NCOL);
  junk = malloc(sizeof *junk * NCOL);
  reunit_work = malloc(sizeof *reunit_work * 6 * NCOL);
  reunit_Rwork = malloc(sizeof *reunit_Rwork * 5 * NCOL);

  // Do whatever is needed to get lattice
  startlat_p = reload_lattice(startflag, startfile);

  return 0;
}
// -----------------------------------------------------------------
