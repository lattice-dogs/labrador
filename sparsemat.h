#ifndef SPARSEMAT_H
#define SPARSEMAT_H

#include <stddef.h>
#include <stdint.h>
#include "polx.h"

typedef struct {
  size_t len;
  size_t *rows;
  size_t *cols;
  polx *coeffs;
} sparsemat;

void sparsemat_polx_mul(sparsemat *r, const polx *c, const sparsemat *a);
void sparsemat_polx_mul_add(sparsemat *r, const polx *c, const sparsemat *a);
void sparsemat_scale_add(sparsemat *r, const sparsemat *a, int64_t s);

#endif
