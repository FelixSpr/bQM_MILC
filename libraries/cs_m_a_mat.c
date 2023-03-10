// -----------------------------------------------------------------
// Add result of complex scalar multiplication on matrix
// c <-- c + s * b
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/bQM.h"

void c_scalar_mult_sum_mat(matrix *b, complex *s, matrix *c) {
  register int i, j;
  for (i = 0; i < NCOL; i++) {
    for (j = 0; j < NCOL; j++) {
      c->e[i][j].real += b->e[i][j].real * s->real - b->e[i][j].imag * s->imag;
      c->e[i][j].imag += b->e[i][j].imag * s->real + b->e[i][j].real * s->imag;
    }
  }
}
// -----------------------------------------------------------------
