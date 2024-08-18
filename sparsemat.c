#include <stddef.h>
#include <stdint.h>
#include "polx.h"
#include "sparsemat.h"

void sparsemat_polx_mul(sparsemat *r, const polx *c, const sparsemat *a) {
  size_t i;

  r->len = a->len;
  for(i=0;i<a->len;i++) {
    r->rows[i] = a->rows[i];
    r->cols[i] = a->cols[i];
    polx_mul(&r->coeffs[i],c,&a->coeffs[i]);
  }
}

void sparsemat_polx_mul_add(sparsemat *r, const polx *c, const sparsemat *a) {
  size_t i,j;

  for(i=0;i<a->len;i++) {
    for(j=0;j<r->len;j++) {
      if(r->rows[j] == a->rows[i] && r->cols[j] == a->cols[i]) {
        polx_mul_add(&r->coeffs[j],c,&a->coeffs[i]);
        break;
      }
    }
    if(j == r->len) {
      r->rows[j] = a->rows[i];
      r->cols[j] = a->cols[i];
      polx_mul(&r->coeffs[j],c,&a->coeffs[i]);
      r->len += 1;
    }
  }
}

void sparsemat_scale_add(sparsemat *r, const sparsemat *a, int64_t s) {
  size_t i,j;

  for(i=0;i<a->len;i++) {
    for(j=0;j<r->len;j++) {
      if(r->rows[j] == a->rows[i] && r->cols[j] == a->cols[i]) {
        polx_scale_add(&r->coeffs[j],&a->coeffs[i],s);
        break;
      }
    }
    if(j == r->len) {
      r->rows[j] = a->rows[i];
      r->cols[j] = a->cols[i];
      polx_scale(&r->coeffs[j],&a->coeffs[i],s);
      r->len += 1;
    }
  }
}
