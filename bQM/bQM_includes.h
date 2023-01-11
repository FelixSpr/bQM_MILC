// -----------------------------------------------------------------
// Include files for supersymmetric evolution
#include "../include/config.h"  // Keep this first
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>             // For setup.c, gauge_info.c
#include "../include/complex.h"
#include "../include/bQM.h"
#include "../include/macros.h"
#include "lattice.h"
#include "../include/comdefs.h"
#include "../include/io_lat.h"
#include "../include/generic.h"
#include "../include/dirs.h"
#include "../include/field_alloc.h"
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Prototypes for functions in high level code
int setup();
int readin(int prompt);
void update();
void update_h(Real eps);
void update_u(Real eps);
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Initialization and set up
void setup_lambda();

// Gaussian random momentum matrices and pseudofermions
void ranmom();

// Polyakov loop observables
complex ploop();
complex ploop_eig();
// Use LAPACK to diagonalize Polyakov loop
// http://www.physics.orst.edu/~rubin/nacphy/lapack/routines/zgeev.html
// First two arguments turn off eigenvector computations
// Third and fifth arguments are the dimensions of the matrix
// Fourth argument is that matrix, which will be overwritten
// Sixth argument holds the computed eigenvalues
// Seventh and ninth arguments are eigenvectors
// Eighth and tenth arguments are the dimensions of the eigenvectors
// Eleventh argument is real workspace, of size given by the twelfth argument
// Thirteenth argument is real workspace, of size given by the third argument
// Final argument reports success or information about failure
void zgeev_(char *doL, char *doR, int *N1, double *store, int *N2, double *eigs,
            double *dumL, int *NL, double *dumR, int *NR,
            double *work, int *Nwork, double *Rwork, int *stat);


// Scalar eigenvalues: averages, extrema and width
// #define SCALAR_EIG_DIST to print all eigenvalues in serial
void scalar_eig(double *ave_eigs, double *eig_widths,
                double *min_eigs, double *max_eigs);

// Action routines
double action();
double bosonic_action();

// Force routines
double bosonic_force(Real eps);

// Compute average Tr[X[i] X[i]] / N_c
double scalar_trace(double *Xtr, double *Xwidth);

// Routines in library_util.c that loop over all sites
#ifdef HMC_ALGORITHM
void copy_bosons(int sign);
#endif
void shiftmat(matrix *dat, matrix *temp, int dir);
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// More measurements
// Use LAPACK for the scalar eigenvalues
// http://www.physics.orst.edu/~rubin/nacphy/lapack/routines/zheev.html
// First argument turns on eigenvector computations
// Second argument chooses between storing upper or lower triangle
// Third and fifth arguments are the dimensions of the matrix
// Fourth argument is that matrix, overwritten by the eigenvectors
// Sixth argument holds the computed eigenvalues
// Seventh argument is complex workspace of size given by the eighth argument
// Ninth argument is real workspace of size 3 * NCOL - 2
// Final argument reports success or information about failure
void zheev_(char *doV, char *uplo, int *N1, double *store, int *N2,
            double *eigs, double *work, int *Nwork, double *Rwork, int *stat);
// -----------------------------------------------------------------
