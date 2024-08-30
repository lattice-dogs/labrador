#include <stdint.h>
#include "data.h"
#include "polx.h"
#include "poly.h"
#include "polz.h"

static int64_t cmodq(int64_t a) {
  int64_t t;
  const int64_t mask = ((int64_t)1 << LOGQ) - 1;
  const int64_t q = ((int64_t)1 << LOGQ) - QOFF;

  t = a >> LOGQ;
  a &= mask;
  a += t*QOFF;
  t = q/2 - a;
  a -= (t >> 63)&q;
  return a;
}

void polx_print(const polx *a) {
  polz t[1];
  polz_frompolx(t,a);
  polz_center(t);
  polz_print(t);
}

void polx_getcoeff(zz *r, const polx *a, int k) {
  polz z;

  polz_frompolx(&z,a);
  polz_getcoeff(r,&z,k);
}

void polxvec_setzero(polx *r, size_t len) {
  polyvec_setzero(&r->vec[0],len*K);
}

void polxvec_copy(polx *r, const polx *a, size_t len) {
  polyvec_copy(&r->vec[0],&a->vec[0],len*K);
}

void polxvec_fromint64vec(polx *r, size_t len, size_t deg, const int64_t a[len*deg*N]) {
  size_t i;
  polz t[deg];

  for(i=0;i<len;i++) {
    polzvec_fromint64vec(t,1,deg,&a[i*deg*N]);
    polzvec_topolxvec(&r[i*deg],t,deg);
  }
}

int polx_iszero(const polx *a) {
  polz t;

  polz_frompolx(&t,a);
  polz_center(&t);
  return polz_iszero(&t);
}

int polxvec_iszero(const polx *a, size_t len) {
  polz t[len];

  polzvec_frompolxvec(t,a,len);
  polzvec_center(t,len);
  return polzvec_iszero(t,len);
}

int polx_iszero_constcoeff(const polx *a) {
  polz t;

  polz_frompolx(&t,a);
  polz_center(&t);
  return polz_iszero_constcoeff(&t);
}

void polx_monomial(polx *r, int64_t v, int k) {
  size_t i;

  v = cmodq(v);
  for(i=0;i<K;i++)
    poly_monomial_ntt(&r->vec[i],v % primes[i].p,k,&primes[i]);
}

void polxvec_ternary(polx *r, size_t len, const uint8_t seed[16], uint64_t nonce) {
  poly t[32];

  while(len >= 32) {
    polyvec_ternary(t,32,seed,nonce);
    polxvec_frompolyvec(r,t,32);
    nonce += (uint64_t)1 << 32;
    r += 32;
    len -= 32;
  }

  polyvec_ternary(t,len,seed,nonce);
  polxvec_frompolyvec(r,t,len);
}

void polxvec_quarternary(polx *r, size_t len, const uint8_t seed[16], uint64_t nonce) {
  poly t[32];

  while(len >= 32) {
    polyvec_quarternary(t,32,seed,nonce);
    polxvec_frompolyvec(r,t,32);
    nonce += (uint64_t)1 << 32;
    r += 32;
    len -= 32;
  }

  polyvec_quarternary(t,len,seed,nonce);
  polxvec_frompolyvec(r,t,len);
}

void polxvec_challenge(polx *r, size_t len, const uint8_t seed[16], uint64_t nonce) {
  poly t[10];

  while(len >= 10) {
    polyvec_challenge(t,10,seed,nonce);
    polxvec_frompolyvec(r,t,10);
    nonce += (uint64_t)1 << 32;
    r += 10;
    len -= 10;
  }

  polyvec_challenge(t,len,seed,nonce);
  polxvec_frompolyvec(r,t,len);
}

void polxvec_almostuniform(polx *r, size_t len, const uint8_t seed[16], uint64_t nonce) {
  polz t[32];

  while(len >= 32) {
    polzvec_almostuniform(t,32,seed,nonce);
    polzvec_topolxvec(r,t,32);
    nonce += (uint64_t)1 << 32;
    r += 32;
    len -= 32;
  }

  polzvec_almostuniform(t,len,seed,nonce);
  polzvec_topolxvec(r,t,len);
}

void polx_frompoly(polx *r, const poly *a) {
  size_t i;

  for(i=0;i<K;i++)
    poly_scale(&r->vec[i],a,primes[i].s,&primes[i]);

  polx_ntt(r,r);
}

void polxvec_frompolyvec(polx *r, const poly *a, size_t len) {
  size_t i;

  for(i=0;i<K;i++)
    polyvec_scale_widening(&r->vec[i],a,len,K,primes[i].s,&primes[i]);

  polxvec_ntt(r,r,len);
}

void polx_refresh(polx *r) {
  polz a;

  polz_frompolx(&a,r);
  polz_center(&a);
  polz_topolx(r,&a);
}

void polxvec_refresh(polx *r, size_t len) {
  size_t i;

  for(i=0;i<len;i++)
    polx_refresh(&r[i]);
}

void polx_reduce(polx *r) {
  size_t i;

  for(i=0;i<K;i++)
    poly_reduce(&r->vec[i],&primes[i]);
}

void polxvec_reduce(polx *r, size_t len) {
  size_t i;

  for(i=0;i<K;i++)
    polyvec_reduce(&r->vec[i],len,K,&primes[i]);
}

void polx_neg(polx *r, const polx *a) {
  polyvec_neg(&r->vec[0],&a->vec[0],K);
}

void polxvec_neg(polx *r, const polx *a, size_t len) {
  polyvec_neg(&r->vec[0],&a->vec[0],K*len);
}

void polx_add(polx *r, const polx *a, const polx *b) {
  polyvec_add(&r->vec[0],&a->vec[0],&b->vec[0],K);
  polx_reduce(r);
}

void polxvec_add(polx *r, const polx *a, const polx *b, size_t len) {
  polyvec_add(&r->vec[0],&a->vec[0],&b->vec[0],K*len);
  polxvec_reduce(r,len);
}

void polx_sub(polx *r,const polx *a,const polx *b) {
  polyvec_sub(&r->vec[0],&a->vec[0],&b->vec[0],K);
  polx_reduce(r);
}

void polxvec_sub(polx *r, const polx *a, const polx *b, size_t len) {
  polyvec_sub(&r->vec[0],&a->vec[0],&b->vec[0],K*len);
  polxvec_reduce(r,len);
}

void polx_ntt(polx *r, const polx *a) {
  size_t i;

  for(i=0;i<K;i++)
    poly_ntt(&r->vec[i],&a->vec[i],&primes[i]);
}

void polxvec_ntt(polx *r, const polx *a, size_t len) {
  size_t i;

  for(i=0;i<K;i++)
    polyvec_ntt(&r->vec[i],&a->vec[i],len,K,&primes[i]);
}

void polx_invntt(polx *r, const polx *a) {
  size_t i;

  for(i=0;i<K;i++)
    poly_invntt(&r->vec[i],&a->vec[i],&primes[i]);
}

void polxvec_invntt(polx *r, const polx *a, size_t len) {
  size_t i;

  for(i=0;i<K;i++)
    polyvec_invntt(&r->vec[i],&a->vec[i],len,K,&primes[i]);
}

void polx_mul(polx *r, const polx *a, const polx *b) {
  size_t i;

  for(i=0;i<K;i++)
    poly_pointwise(&r->vec[i],&a->vec[i],&b->vec[i],&primes[i]);
}

void polx_poly_mul(polx *r, const polx *a, const poly *b) {
  polx t;

  polx_frompoly(&t,b);
  polx_mul(r,a,&t);
}

void polxvec_mul(polx *r, const polx *a, const polx *b, size_t len) {
  size_t i;

  for(i=0;i<K;i++)
    polyvec_pointwise(&r->vec[i],&a->vec[i],&b->vec[i],len,K,&primes[i]);
}

void polxvec_polx_mul(polx *r, const polx *a, const polx *b, size_t len) {
  size_t i;

  for(i=0;i<K;i++)
    polyvec_poly_pointwise(&r->vec[i],&a->vec[i],&b->vec[i],len,K,&primes[i]);
}

void polx_mul_add(polx *r, const polx *a, const polx *b) {
  size_t i;

  for(i=0;i<K;i++)
    poly_pointwise_add(&r->vec[i],&a->vec[i],&b->vec[i],&primes[i]);

  polx_reduce(r);
}

void polxvec_mul_add(polx *r, const polx *a, const polx *b, size_t len) {
  size_t i;

  for(i=0;i<K;i++)
    polyvec_pointwise_add(&r->vec[i],&a->vec[i],&b->vec[i],len,K,&primes[i]);

  polxvec_reduce(r,len);
}

void polxvec_polx_mul_add(polx *r, const polx *a, const polx *b, size_t len) {
  size_t i;

  for(i=0;i<K;i++)
    polyvec_poly_pointwise_add(&r->vec[i],&a->vec[i],&b->vec[i],len,K,&primes[i]);

  polxvec_reduce(r,len);
}

void polxvec_sprod(polx *r, const polx *a, const polx *b, size_t len) {
  size_t i;

  if(len == 0) {
    polxvec_setzero(r,1);
    return;
  }

  for(i=0;i<K;i++)
    polyvec_sprod_pointwise(&r->vec[i],&a->vec[i],&b->vec[i],len,K,&primes[i]);
}

void polxvec_sprod_add(polx *r, const polx *a, const polx *b, size_t len) {
  size_t i;

  for(i=0;i<K;i++)
    polyvec_sprod_pointwise_add(&r->vec[i],&a->vec[i],&b->vec[i],len,K,&primes[i]);
}

void polx_scale(polx *r, const polx *a, int64_t s) {
  size_t i;

  s = cmodq(s);
  s <<= 16;
  for(i=0;i<K;i++)
    poly_scale(&r->vec[i],&a->vec[i],s % primes[i].p,&primes[i]);
}

void polx_scale_frompoly(polx *r, const poly *a, int64_t s) {
  size_t i;

  s = cmodq(s);
  for(i=0;i<K;i++)
    poly_scale(&r->vec[i],a,s*primes[i].s % primes[i].p,&primes[i]);

  polx_ntt(r,r);
}

void polxvec_scale(polx *r, const polx *a, size_t len, int64_t s) {
  size_t i;

  s = cmodq(s);
  s <<= 16;
  for(i=0;i<K;i++)
    polyvec_scale(&r->vec[i],&a->vec[i],len,K,s % primes[i].p,&primes[i]);
}

void polxvec_scale_frompolyvec(polx *r, const poly *a, size_t len, int64_t s) {
  size_t i;

  s = cmodq(s);
  for(i=0;i<K;i++)
    polyvec_scale_widening(&r->vec[i],a,len,K,s*primes[i].s % primes[i].p,&primes[i]);

  polxvec_ntt(r,r,len);
}

void polx_scale_add(polx *r, const polx *a, int64_t s) {
  size_t i;

  s = cmodq(s);
  s <<= 16;
  for(i=0;i<K;i++)
    poly_scale_add(&r->vec[i],&a->vec[i],s % primes[i].p,&primes[i]);

  polx_reduce(r);
}

void polxvec_scale_add(polx *r, const polx *a, size_t len, int64_t s) {
  size_t i;

  s = cmodq(s);
  s <<= 16;
  for(i=0;i<K;i++)
    polyvec_scale_add(&r->vec[i],&a->vec[i],len,K,s % primes[i].p,&primes[i]);

  polxvec_reduce(r,len);
}

size_t polxvec_mul_extension(polx *c, const polx *a, const polx *b, size_t len, size_t deg, size_t mult) {
  size_t i,j,k = 0;

  for(i=0;i<mult;i++)
    for(j=0;j<K;j++)
      k = polyvec_pointwise_extension(&c[i].vec[j],&a->vec[j],&b[i].vec[j],len,K*mult,deg,&primes[j]);

  return k;
}

size_t polxvec_collaps_add_extension(polx *c, const polx *a, const polx *b, size_t len, size_t deg, size_t mult) {
  size_t i,j,k = 0;

  for(i=0;i<mult;i++)
    for(j=0;j<K;j++)
      k = polyvec_collaps_add_extension(&c[i].vec[j],&a[i].vec[j],&b->vec[j],len,K*mult,deg,&primes[j]);

  return k;
}

void polxvec_decompose(poly *r, const polx *a, size_t len, size_t t, size_t d) {
  size_t i;
  const size_t stride = len;
  polz b[16];

  while(len >= 16) {
    polzvec_frompolxvec(b,a,16);
    polzvec_center(b,16);
    for(i=0;i<16;i++)
      polz_decompose(&r[i],&b[i],stride,t,d);
    r += 16;
    a += 16;
    len -= 16;
  }

  polzvec_frompolxvec(b,a,len);
  polzvec_center(b,len);
  for(i=0;i<len;i++)
    polz_decompose(&r[i],&b[i],stride,t,d);
}

void polxvec_reconstruct(polx *r, const poly *a, size_t len, size_t t, size_t d) {
  size_t i;
  const size_t stride = len;
  polz b[16];

  while(len >= 16) {
    for(i=0;i<16;i++)
      polz_reconstruct(&b[i],&a[i],stride,t,d);
    polzvec_topolxvec(r,b,16);
    r += 16;
    a += 16;
    len -= 16;
  }

  for(i=0;i<len;i++)
    polz_reconstruct(&b[i],&a[i],stride,t,d);
  polzvec_topolxvec(r,b,len);
}

void polx_sigmam1(polx *r, const polx *a) {
  polyvec_sigmam1_ntt(&r->vec[0],&a->vec[0],K);
}

void polxvec_sigmam1(polx *r, const polx *a, size_t len) {
  polyvec_sigmam1_ntt(&r->vec[0],&a->vec[0],K*len);
}

/* FIXME: NTT rep! */
void polx_sigma5(polx *r, const polx *a) {
  polyvec_sigma5(&r->vec[0],&a->vec[0],K);
}

void polxvec_sigma5(polx *r, const polx *a, size_t len) {
  polyvec_sigma5(&r->vec[0],&a->vec[0],K*len);
}

void polx_sigma5inv(polx *r, const polx *a) {
  polyvec_sigma5inv(&r->vec[0],&a->vec[0],K);
}

void polxvec_sigma5inv(polx *r, const polx *a, size_t len) {
  polyvec_sigma5inv(&r->vec[0],&a->vec[0],K*len);
}

void polx_flip(polx *r, const polx *a) {
  size_t i;

  for(i=0;i<K;i++)
    poly_flip_ntt(&r->vec[i],&a->vec[i],&primes[i]);
}

void polxvec_flip(polx *r, const polx *a, size_t len) {
  size_t i;

  for(i=0;i<K;i++)
    polyvec_flip_ntt(&r->vec[i],&a->vec[i],len,K,&primes[i]);
}

