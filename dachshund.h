#ifndef DACHSHUND_H
#define DACHSHUND_H

#include <stdint.h>
#include "poly.h"
#include "polx.h"
#include "polz.h"
#include "labrador.h"
#include "chihuahua.h"

typedef struct {
  size_t r;
  size_t *n;
  size_t fu;
  size_t bu;
  size_t kappa;
  size_t kappa1;
  polz *u;
  polx *alpha;
} commitment;

typedef struct {
  size_t r;
  size_t *n;
  size_t k;
  sparsecnst *cnst;
  uint64_t *betasq;
  uint8_t h[16];
  polx *u0;
  polx *alpha;
} smplstmnt;

/* allocate statement for r witness vectors of ranks n[r] and l2-norm bounds betasq[r] (0 for 0/1),
 * consisting of k linear constraints (rows of constraint matrix) */
#define init_smplstmnt_raw NAMESPACE(init_smplstmnt_raw)
__attribute__((visibility("default")))
void init_smplstmnt_raw(smplstmnt *st, size_t r, const size_t n[r], const uint64_t betasq[r], size_t k);
/* set j-th linear constraint (j <= k) where the nz sub-vectors of indices idx[nz] (0 <= idx[nz] < r) are non-zero;
 * phi is the concatenation of the non-zero sub vectors (of length \sum_{i=0}^{nz-1} n[i]) */
#define set_smplstmnt_lincnst_raw NAMESPACE(set_smplstmnt_lincnst_raw)
__attribute__((visibility("default")))
int set_smplstmnt_lincnst_raw(smplstmnt *st, size_t i, size_t nz, const size_t idx[nz],
                              const size_t n[nz], size_t deg, int64_t *phi, int64_t *b);
void print_smplstmnt_pp(const smplstmnt *st);
#define free_smplstmnt NAMESPACE(free_smplstmnt)
__attribute__((visibility("default")))
void free_smplstmnt(smplstmnt *st);
#define free_commitment NAMESPACE(free_commitment)
__attribute__((visibility("default")))
void free_commitment(commitment *com);

void simple_prove(statement *ost, witness *owt, proof *pi, commitment *com,
                  const smplstmnt *ist, const witness *iwt, int tail);
int simple_reduce(statement *ost, const proof *pi, const commitment *com, const smplstmnt *ist);
#define simple_verify NAMESPACE(simple_verify)
__attribute__((visibility("default")))
int simple_verify(const smplstmnt *st, const witness *wt);

#endif
