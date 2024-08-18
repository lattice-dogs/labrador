#ifndef PACK_H
#define PACK_H

#include <stdint.h>
#include <stddef.h>
#include "poly.h"
#include "polx.h"
#include "polz.h"
#include "labrador.h"
#include "chihuahua.h"
#include "dachshund.h"
#include "greyhound.h"

typedef struct {
  size_t l;
  double size;
  proof *pi[16];
  witness owt;
} composite;

#define free_composite NAMESPACE(free_composite)
__attribute__((visibility("default")))
void free_composite(composite *proof);

double composite_prove_principle(composite *proof, const prncplstmnt *st, const witness *wt);
#define composite_prove_simple NAMESPACE(composite_prove_simple)
__attribute__((visibility("default")))
double composite_prove_simple(composite *proof, commitment *com, const smplstmnt *st, const witness *wt);
double composite_prove_polcom(composite *proof, polcomprf *ppi, polcomctx *ctx, uint32_t x, uint32_t y);

int composite_verify_principle(const composite *proof, const prncplstmnt *st);
#define composite_verify_simple NAMESPACE(composite_verify_simple)
__attribute__((visibility("default")))
int composite_verify_simple(const composite *proof, const commitment *com, const smplstmnt *st);
int composite_verify_polcom(const composite *proof, const polcomprf *ppi);

#endif
