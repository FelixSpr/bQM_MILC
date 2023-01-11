// -----------------------------------------------------------------
// Define global scalars and fields in the lattice
#ifndef _LATTICE_H
#define _LATTICE_H

#include "defines.h"
#include "../include/macros.h"    // For MAXFILENAME
#include "../include/io_lat.h"    // For gauge_file
#include "../include/bQM.h"
#include "../include/random.h"    // For double_prn
#include "../include/dirs.h"      // For NDIMS
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// The lattice is an array of this site struct
typedef struct {
  short t;            // Coordinates of this site
  char parity;        // Is it even or odd?
  int index;          // Index in the array

#ifdef SITERAND
  // The state information for a random number generator
  double_prn site_prn;
#endif

  matrix link;          // Gauge link (in group)
  matrix X[NSCALAR];    // Scalars (in algebra)

#ifdef HMC_ALGORITHM
  matrix old_link, old_X[NSCALAR];      // For accept/reject
#endif

  // All momenta should be anti-hermitian matrices
  // Since scalars are anti-hermitian but stored as full matrices,
  // simplest to treat their momenta the same way
  anti_hermitmat mom;
  matrix f_U, mom_X[NSCALAR], f_X[NSCALAR];
} site;
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Definition of global variables
#ifdef CONTROL
#define EXTERN
#else
#define EXTERN extern
#endif

EXTERN int nt;              // Lattice length
EXTERN int iseed;           // Random number seed
EXTERN int warms, trajecs, propinterval;
EXTERN Real traj_length;

// SU(N) generators
EXTERN matrix Lambda[NADJ];
EXTERN Real beta, omega;
EXTERN double_complex linktrsum;
EXTERN u_int32type nersc_checksum;
EXTERN char startfile[MAXFILENAME], savefile[MAXFILENAME];
EXTERN int startflag;     // Beginning lattice: CONTINUE, RELOAD, FRESH
EXTERN int fixflag;       // Gauge fixing: COULOMB_GAUGE_FIX, NO_GAUGE_FIX
EXTERN int saveflag;      // 1 if we will save the lattice;

// Some of these global variables are node dependent
// They are set in "make_lattice()"
EXTERN int sites_on_node;       // Number of sites on this node
EXTERN int even_sites_on_node;  // Number of even sites on this node
EXTERN int odd_sites_on_node;   // Number of odd sites on this node
EXTERN int number_of_nodes;     // Number of nodes in use
EXTERN int this_node;           // Node number of this node

// Stuff for HMC
EXTERN int nsteps;
EXTERN Real bnorm, max_bf;

// Each node maintains a structure with the pseudorandom number
// generator state
EXTERN double_prn node_prn;

// For convenience in calculating action and force
EXTERN Real one_ov_N;

// Temporary matrices
EXTERN matrix *tempmat, *tempmat2, *temp_X[NSCALAR];

// Allocate some more arrays to be used by LAPACK in scalar eigenvalues
// and reunitarization (in generic directory)
EXTERN double *Rwork, *eigs, *store, *work;
EXTERN double *junk, *left, *right;           // For reunitarization
EXTERN double *reunit_work, *reunit_Rwork;    // For reunitarization

EXTERN gauge_file *startlat_p;
EXTERN gauge_file *savelat_p;

// The lattice is a single global variable
// (actually this is the part of the lattice on this node)
EXTERN site *lattice;

// Vectors for addressing
// Generic pointers, for gather routines
// Need at least 2NSCALAR+1=19 to gather scalars for the bosonic force
#define N_POINTERS 20
EXTERN char **gen_pt[N_POINTERS];

#endif // _LATTICE_H
// -----------------------------------------------------------------
