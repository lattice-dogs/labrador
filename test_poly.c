#include <stdio.h>
#include <math.h>
#include "data.h"
#include "randombytes.h"
#include "poly.h"

static void poly_naivemul(poly *c,const poly *a,const poly *b,const pdata *prime) {
  int i,j;
  int32_t r[2*N] = {0};

  for(i=0;i<N;i++)
    for(j=0;j<N;j++)
       r[i+j] += (int32_t)a->vec->c[i]*b->vec->c[j];

  for(i=0;i<N;i++)
    c->vec->c[i] = (r[i] - r[N+i]) % prime->p;
}

int main(void) {
  int i;
  const pdata *prime = &primes[K-1];
  __attribute__((aligned(16)))
  uint8_t seed[16];
  poly a,b,c,d;

  randombytes(seed,16);
  polyvec_uniform(&a,1,prime,seed,0);
  polyvec_ternary(&b,1,seed,1);
  polyvec_challenge(&c,1,seed,2);

  printf("norm uniform: %.2f (expected: %.2f)\n",polyvec_norm(&a,1),prime->p*sqrt(N/12.0));
  printf("norm ternary: %.2f (expected: %.2f)\n",polyvec_norm(&b,1),sqrt(10*N/16));
  printf("norm challenge: %.2f (expected: %.2f)\n",polyvec_norm(&c,1),sqrt(TAU1+4*TAU2));

  poly_naivemul(&d,&a,&b,prime);

/*
  poly_scale(&a,&a,prime->s,prime);
  poly_scale(&b,&b,prime->s,prime);
  poly_ntt(&a,&a,prime);
  poly_ntt(&b,&b,prime);
  poly_pointwise(&c,&a,&b,prime);
  poly_invntt(&c,&c,prime);
  poly_scale(&c,&c,,prime);
  poly_sub(&d,&d,&c);
  poly_reduce(&d,prime);
  for(i=0;i<N;i++)
    if(d.vec->c[i])
      fprintf(stderr,"ERROR in mul: %d %d\n", i, d.vec->c[i]);
*/

  poly_sigmam1(&c,&a);
  poly_sigma(&d,&a,-1);
  poly_sub(&c,&c,&d);
  poly_reduce(&c,prime);
  for(i=0;i<N;i++)
    if(c.vec->c[i])
      fprintf(stderr,"ERROR in poly_sigmam1: %i\n", i);

  poly_sigmam1(&c,&a);
  poly_sigmam1(&c,&c);
  poly_sub(&c,&c,&a);
  poly_reduce(&c,prime);
  for(i=0;i<N;i++)
    if(c.vec->c[i])
      fprintf(stderr,"ERROR in poly_sigmam1: %i\n", i);

  poly_sigma5(&c,&a);
  poly_sigma5inv(&c,&c);
  poly_sub(&c,&c,&a);
  poly_reduce(&c,prime);
  for(i=0;i<N;i++)
    if(c.vec->c[i])
      fprintf(stderr,"ERROR in poly_sigma5/sigma5inv: %i\n", i);

  return 0;
}
