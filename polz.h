#ifndef POLZ_H
#define POLZ_H

#include <stdint.h>
#include <immintrin.h>
#include "data.h"
#include "polx.h"
#include "poly.h"

typedef struct {
  vecn limbs[L];
} polz;

void polz_print(const polz *a);
void polzvec_copy(polz *r, const polz *a, size_t len);
void polz_getcoeff(zz *r, const polz *a, int k);
void polz_setcoeff(polz *r, const zz *a, int k);
void polz_setcoeff_fromint64(polz *r, int64_t a, int k);
void polzvec_fromint64vec(polz *r, size_t len, size_t deg, const int64_t v[len*deg*N]);
int polz_iszero(const polz *a);
int polz_iszero_constcoeff(const polz *a);
int polzvec_iszero(const polz *a, size_t len);

void polzvec_uniform(polz *r, size_t len, const uint8_t seed[16], uint64_t nonce);
void polz_bitpack(uint8_t r[N*QBYTES], const polz *a);
void polzvec_bitpack(uint8_t *r, const polz *a, size_t len);
void polz_bitunpack(polz *r, const uint8_t buf[N*QBYTES]);
void polzvec_almostuniform(polz *r, size_t len, const uint8_t seed[16], uint64_t nonce);
double polzvec_norm(const polz *r, size_t len);

void polz_reduce(polz *r);
void polzvec_reduce(polz *r, size_t len);
void polz_caddq(polz *r);
void polzvec_caddq(polz *r, size_t len);
void polz_center(polz *r);
void polzvec_center(polz *r, size_t len);

void polz_topoly(poly *r, const polz *a);
void polzvec_topolyvec(poly *r, const polz *a, size_t len, size_t stride);
void polz_frompoly(polz *r, const poly *a);
void polzvec_frompolyvec(polz *r, const poly *a, size_t len);
void polz_topolx(polx *r, const polz *a);
void polzvec_topolxvec(polx *r, const polz *a, size_t len);
void polz_frompolx(polz *r, const polx *a);
void polzvec_frompolxvec(polz *r, const polx *a, size_t len);

void polz_add(polz *r, const polz *a, const polz *b);
void polzvec_add(polz *r, const polz *a, const polz *b, size_t len);
void polz_sub(polz *r, const polz *a, const polz *b);
void polzvec_sub(polz *r, const polz *a, const polz *b, size_t len);
void polz_slli(polz *r, const polz *a, int s);
void polzvec_slli(polz *r, const polz *a, size_t len, int s);

void polz_mul(polz *r, const polz *a, const polz *b);
void polz_poly_mul(polz *r, const polz *a, const poly *b);

void polz_split(poly *l, polz *h, const polz *a, size_t d);
void polzvec_split(poly *l, polz *h, const polz *a, size_t len, size_t d);
void polz_decompose(poly *l, const polz *a, size_t stride, size_t t, size_t d);
void polzvec_decompose(poly *r, const polz *a, size_t len, size_t t, size_t d);
void polz_decompose_topolx(polx *r, const polz *a, size_t stride, size_t t, size_t d);
void polzvec_decompose_topolxvec(polx *r, const polz *a, size_t len, size_t stride, size_t t, size_t d);
void polz_reconstruct(polz *r, const poly *a, size_t stride, size_t t, size_t d);
void polzvec_reconstruct(polz *r, const poly *a, size_t len, size_t t, size_t d);

void polz_sigmam1(polz *r, const polz *a);
void polzvec_sigmam1(polz *r, const polz *a, size_t len);

#endif
