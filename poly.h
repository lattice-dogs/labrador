#ifndef POLY_H
#define POLY_H

#include <stdint.h>
#include <immintrin.h>
#include <complex.h>
#include "data.h"

typedef union {
  __m512i v[N/32];
  int16_t c[N];
} vecn;

typedef struct {
  vecn vec[1];
} poly;

__attribute__((const))
size_t extlen(size_t len, size_t deg);

void polyvec_setzero(poly *r, size_t len);
int polyvec_isbinary(const poly *r, size_t len);
void polyvec_fromint64vec(poly *r, size_t len, size_t deg, const int64_t a[len*deg*N]);
void polyvec_copy(poly *r, const poly *a, size_t len);
void poly_binary_fromuint64(poly *r, uint64_t a);
void poly_monomial_ntt(poly *r, int16_t v, int k, const pdata *prime);
void polyvec_uniform(poly *r, size_t len, const pdata *prime, const uint8_t seed[16], uint64_t nonce);
void polyvec_ternary(poly *r, size_t len, const uint8_t seed[16], uint64_t nonce);
void polyvec_quarternary(poly *r, size_t len, const uint8_t seed[16], uint64_t nonce);
void polyvec_challenge(poly *c, size_t len, const uint8_t seed[16], uint64_t nonce);

int64_t polyvec_sprodz_ref(const poly *a, const poly *b, size_t len);
int64_t polyvec_sprodz(const poly *a, const poly *b, size_t len);
double polyvec_norm(const poly *a, size_t len);

void poly_reduce(poly *r, const pdata *prime);
void polyvec_reduce(poly *r, size_t len, size_t stride, const pdata *prime);
void poly_center(poly *r, const pdata *prime);
void polyvec_center(poly *r, size_t len, size_t stride, const pdata *prime);
void poly_csubp(poly *r, const pdata *prime);
void polyvec_csubp(poly *r, size_t len, size_t stride, const pdata *prime);
void poly_caddp(poly *r, const pdata *prime);
void polyvec_caddp(poly *r, size_t len, size_t stride, const pdata *prime);
void poly_quot_add(poly *r, const poly *a, const pdata *prime);
void polyvec_quot_add(poly *r, const poly *a, size_t len, size_t stride, const pdata *prime);

void poly_neg(poly *r, const poly *a);
void polyvec_neg(poly *r, const poly *a, size_t len);
void poly_add(poly *r, const poly *a, const poly *b);
void polyvec_add(poly *r, const poly *a, const poly *b, size_t len);
void poly_sub(poly *r, const poly *a, const poly *b);
void polyvec_sub(poly *r, const poly *a, const poly *b, size_t len);

void poly_ntt_ref(poly *a, const pdata *prime);
__attribute__((visibility("hidden")))
void poly_ntt(poly *r, const poly *a, const pdata *prime);
void polyvec_ntt(poly *r, const poly *a, size_t len, size_t stride, const pdata *prime);
__attribute__((visibility("hidden")))
void poly_invntt(poly *r, const poly *a, const pdata *prime);
void polyvec_invntt(poly *r, const poly *a, size_t len, size_t stride, const pdata *prime);

void poly_pointwise(poly *r, const poly *a, const poly *b, const pdata *prime);
void polyvec_pointwise(poly *r, const poly *a, const poly *b, size_t len, size_t stride, const pdata *prime);
void polyvec_poly_pointwise(poly *r, const poly *a, const poly *b, size_t len, size_t stride, const pdata *prime);
void poly_pointwise_add(poly *r, const poly *a, const poly *b, const pdata *prime);
void polyvec_pointwise_add(poly *r, const poly *a, const poly *b, size_t len, size_t stride, const pdata *prime);
void polyvec_poly_pointwise_add(poly *r, const poly *a, const poly *b, size_t len, size_t stride, const pdata *prime);
void polyvec_sprod_pointwise(poly *r, const poly *a, const poly *b, size_t len, size_t stride, const pdata *prime);
void polyvec_sprod_pointwise_add(poly *r, const poly *a, const poly *b, size_t len, size_t stride, const pdata *prime);

size_t polyvec_pointwise_extension(poly *c, const poly *a, const poly *b, size_t len, size_t stride,
                                   size_t deg, const pdata *prime);
size_t polyvec_collaps_add_extension(poly *c, const poly *a, const poly *b, size_t len, size_t stride, size_t deg,
                                     const pdata *prime);

void poly_scale(poly *r, const poly *a, int16_t s, const pdata *prime);
void polyvec_scale(poly *r, const poly *a, size_t len, size_t stride, int16_t s, const pdata *prime);
void polyvec_scale_widening(poly *r, const poly *a, size_t len, size_t stride, int16_t s, const pdata *prime);
void poly_scale_add(poly *r, const poly *a, int16_t s, const pdata *prime);
void polyvec_scale_add(poly *r, const poly *a, size_t len, size_t stride, int16_t s, const pdata *prime);

void poly_fft(double complex r[N/2], const poly *a);
void poly_invfft(poly *r, double complex a[N/2]);
double poly_opnorm(const poly *a);

void poly_sigmam1(poly *r, const poly *a);
void polyvec_sigmam1(poly *r, const poly *a, size_t len);
void poly_sigmam1_ntt(poly *r, const poly *a);
void polyvec_sigmam1_ntt(poly *r, const poly *a, size_t len);
void poly_sigma(poly *r, const poly *a, int k);
void poly_sigma5(poly *r, const poly *a);
void polyvec_sigma5(poly *r, const poly *a, size_t len);
void poly_sigma5inv(poly *r, const poly *a);
void polyvec_sigma5inv(poly *r, const poly *a, size_t len);
void poly_flip(poly *r, const poly *a);
void polyvec_flip(poly *r, const poly *a, size_t len);
void poly_flip_ntt(poly *r, const poly *a, const pdata *prime);
void polyvec_flip_ntt(poly *r, const poly *a, size_t len, size_t stride, const pdata *prime);

#endif
