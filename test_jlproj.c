#include <stdint.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "cpucycles.h"
#include "data.h"
#include "randombytes.h"
#include "jlproj.h"
#include "poly.h"
#include "polz.h"

int main(void) {
  size_t i;
  int64_t c;
  unsigned long long t[21], overhead;
  __attribute__((aligned(16)))
  uint8_t seed[16];
  uint64_t nonce = 0;
  __attribute__((aligned(64)))
  uint8_t mat[256*N/8];
  __attribute__((aligned(64)))
  uint8_t buf[256*QBYTES];
  __attribute__((aligned(64)))
  int32_t p[256];
  poly s;
  polx r;
  zz x,y;

  overhead = cpucycles_overhead();

  randombytes(seed,16);
  randombytes(mat,sizeof(mat));
  randombytes(buf,sizeof(buf));
  polyvec_uniform(&s,1,&primes[0],seed,nonce++);

  memset(p,0,sizeof(p));
  poly_jlproj_add(p,&s,mat);
  printf("Projected norm: %.2f; expected: %.2f\n",sqrt(jlproj_normsq(p)),16*polyvec_norm(&s,1));

  polxvec_jlproj_collapsmat(&r,mat,1,buf);
  polx_poly_mul(&r,&r,&s);
  polx_getcoeff(&x,&r,0);
  if(!zz_less_than(&x,&modulus.q))
    zz_sub(&x,&x,&modulus.q);

  c = jlproj_collapsproj(p,buf);
  zz_fromint64(&y,c);
  if(!zz_less_than(&y,&modulus.q))
    zz_sub(&y,&y,&modulus.q);

  if(!zz_equal(&x,&y)) {
    fprintf(stderr,"ERROR: Constant coeff doesn't match\n");
    return 1;
  }

  for(i=0;i<21;i++) {
    t[i] = cpucycles();
    poly_jlproj_add(p,&s,mat);
  }
  for(i=0;i<20;i++)
    printf("poly_jlproj:  %2lu: %8lld\n", i, t[i+1] - t[i] - overhead);

  for(i=0;i<21;i++) {
    t[i] = cpucycles();
    polxvec_jlproj_collapsmat(&r,mat,1,buf);
  }
  for(i=0;i<20;i++)
    printf("jlproj_collapsmat:  %2lu: %8lld\n", i, t[i+1] - t[i] - overhead);

  return 0;
}
