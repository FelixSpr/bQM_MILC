// -----------------------------------------------------------------
// Structure for passing simulation parameters to each node
#ifndef _PARAMS_H
#define _PARAMS_H
#include "../include/macros.h"  // For MAXFILENAME

typedef struct {
  int stopflag;           // 1 if it is time to stop

  // Initialization parameters
  int nt;                 // Lattice dimensions
  int iseed;              // For random numbers

  int warms;              // The number of warmup trajectories
  int trajecs;            // The number of real trajectories
  Real traj_length;       // The length of each trajectory
  int nsteps;             // Steps per trajectory
  int propinterval;       // Number of trajectories between measurements
  int startflag;          // What to do for beginning lattice
  int fixflag;            // Whether to gauge fix to Coulomb gauge
  int saveflag;           // What to do with lattice at end
  Real beta;              // Gauge coupling
  Real omega;             // Quadratic regulator
  char startfile[MAXFILENAME], savefile[MAXFILENAME];
} params;
#endif
// -----------------------------------------------------------------
