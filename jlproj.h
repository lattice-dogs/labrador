#ifndef JLPROJ_H
#define JLPROJ_H

#include <stdint.h>
#include "data.h"
#include "poly.h"
#include "polz.h"

// projection coefficients have stddev sqrt(normsq); must stay below 2^31; hence bound of 2^31/8 = 2^28
#define JLMAXNORM MIN((((uint64_t)1 << LOGQ) - QOFF)/125,(uint64_t)1 << 28)
#define JLMAXNORMSQ (JLMAXNORM*JLMAXNORM)

void poly_jlproj_add(int32_t r[256], const poly *p, const uint8_t mat[256*N/8]);
void polyvec_jlproj_add(int32_t r[256], const poly *p, size_t len, const uint8_t *mat);
void polxvec_jlproj_collapsmat(polx *r, const uint8_t *mat, size_t len, const uint8_t buf[256*QBYTES]);
int64_t jlproj_collapsproj(const int32_t p[256], const uint8_t buf[256*QBYTES]);
uint64_t jlproj_normsq(const int32_t p[256]);

#endif
