// -----------------------------------------------------------------
// Routines for filling momenta with gaussian random numbers
#include "bQM_includes.h"
// -----------------------------------------------------------------



// -----------------------------------------------------------------
// Construct gaussian random momentum matrices
// All need to be anti-hermitian
void ranmom() {
  register int i, j;
  register site *s;
  anti_hermitmat tah;

  FORALLSITES(i, s) {
#ifdef SITERAND
    random_anti_hermitian(&(s->mom), &(s->site_prn));
#else
    random_anti_hermitian(&(s->mom), &(s->node_prn));
#endif

    for (j = 0; j < NSCALAR; j++) {
#ifdef SITERAND
      random_anti_hermitian(&tah, &(s->site_prn));
#else
      random_anti_hermitian(&tah, &(s->node_prn));
#endif
      uncompress_anti_hermitian(&tah, &(s->mom_X[j]));
    }
  }
}
// -----------------------------------------------------------------
