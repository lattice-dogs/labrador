#include <stdint.h>
#include <assert.h>
#include <stdio.h>
#include "randombytes.h"
#include "cpucycles.h"
#include "data.h"
#include "poly.h"

const pdata *prime = &primes[K-1];
static int16_t zeta[N];
static int16_t zetapow[N][N];

static int16_t pow_simple(int16_t a, int e) {
  int16_t r;
  if(e == 0) return 1;
  else if(e == 1) return a;
  else if(e < 0) return pow_simple(a,-e*(prime->p-2));

  r = pow_simple(a,e/2);
  r = (int32_t)r*r % prime->p;
  if(e&1) r = (int32_t)r*a % prime->p;

  return r;
}

int main(void) {
  int i,j;
  unsigned long long t[21], overhead;
  __attribute__((aligned(16)))
  uint8_t seed[16];
  poly a,b;
  int32_t out[N];
  const int16_t s = (int32_t)pow_simple(prime->mont,-1) * prime->s % prime->p;
  const int16_t u = 1024;  // mont*64^-1

  overhead = cpucycles_overhead();
//  for(i=0;i<1000;i++)
//    poly_ntt(&a,prime);

  for(i=0;i<21;i++) {
    t[i] = cpucycles();
    poly_ntt(&a,&a,prime);
  }
  for(i=0;i<20;i++)
    printf("ntt: %2d: %llu\n", i, t[i+1] - t[i] - overhead);
  for(i=0;i<21;i++) {
    t[i] = cpucycles();
    poly_invntt(&a,&a,prime);
  }
  for(i=0;i<20;i++)
    printf("invntt: %2d: %llu\n", i, t[i+1] - t[i] - overhead);

  for(i=0;i<N;i++)
    a.vec->c[i] = 0;
  a.vec->c[1] = 1;
  poly_ntt(&a,&a,prime);
  poly_scale(&a,&a,s,prime);
  for(i=0;i<N;i++) {
    zeta[i] = a.vec->c[i] % prime->p;
    assert((pow_simple(a.vec->c[i],N) + 1) % prime->p == 0);
    for(j=0;j<i;j++)
      assert((zeta[j] - zeta[i]) % prime->p);
  }

  for(i=0;i<N;i++)
    for(j=0;j<N;j++)
      zetapow[i][j] = pow_simple(zeta[i],j);

  for(i=0;i<N;i++) {
    for(j=0;j<N;j++)
      a.vec->c[j] = 0;
    a.vec->c[i] = 1;
    poly_ntt(&a,&a,prime);
    poly_scale(&a,&a,s,prime);
    for(j=0;j<N;j++) {
      a.vec->c[j] %= prime->p;
      assert((a.vec->c[j] - zetapow[j][i]) % prime->p == 0);
    }
  }

  randombytes(seed,16);
  polyvec_uniform(&a,1,prime,seed,0);
  for(i=0;i<N;i++) {
    out[i] = 0;
    for(j=0;j<N;j++)
      out[i] += (int32_t)a.vec->c[j]*zetapow[i][j] % prime->p;
    out[i] %= prime->p;
  }
  poly_ntt(&a,&a,prime);
  poly_scale(&a,&a,s,prime);
  for(i=0;i<N;i++)
    assert((a.vec->c[i] - out[i]) % prime->p == 0);

  polyvec_uniform(&a,1,prime,seed,1);
  poly_ntt(&b,&a,prime);
  poly_scale(&b,&b,s,prime);
  poly_invntt(&b,&b,prime);
  poly_scale(&b,&b,u,prime);
  for(i=0;i<N;i++)
    assert((a.vec->c[i] - b.vec->c[i]) % prime->p == 0);

  for(i=0;i<N;i++)
    a.vec->c[i] = i;
  poly_invntt(&a,&a,prime);
  poly_scale(&a,&a,u,prime);
  poly_sigma5(&a,&a);
  poly_ntt(&a,&a,prime);
  poly_scale(&a,&a,s,prime);
  for(i=0;i<N;i++)
    printf("%d\n",a.vec->c[i]);

  return 0;
}
