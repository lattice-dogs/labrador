#include <stddef.h>
#include <sys/types.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <gmp.h>
#include "data.h"
#include "randombytes.h"
#include "polz.h"
#include "polx.h"
#include "cpucycles.h"

static void polz_kroneckermul(polz *r, const polz *a, const polz *b) {
  ssize_t i,j;
  mpz_t az,bz,qz,tz,uz;
  const mp_bitcnt_t lbeta = 2*LOGQ+6;

  mpz_inits(az,bz,qz,tz,uz,NULL);

  for(i=L-1;i>=0;i--) {
    mpz_mul_2exp(qz,qz,14);
    mpz_add_ui(qz,qz,modulus.q.limbs[i]);
  }

  for(i=N-1;i>=0;i--) {
    for(j=L-1;j>=0;j--) {
      mpz_add_ui(az,az,a->limbs[j].c[i]);
      mpz_add_ui(bz,bz,b->limbs[j].c[i]);
      if(j) {
        mpz_mul_2exp(az,az,14);
        mpz_mul_2exp(bz,bz,14);
      }
    }
    if(i) {
      mpz_mul_2exp(az,az,lbeta-14*(L-1));
      mpz_mul_2exp(bz,bz,lbeta-14*(L-1));
    }
  }

  mpz_mul(az,az,bz);
  mpz_fdiv_q_2exp(bz,az,N*lbeta);

  for(i=0;i<N;i++) {
    mpz_fdiv_r_2exp(tz,az,lbeta);
    mpz_fdiv_q_2exp(az,az,lbeta);
    mpz_fdiv_r_2exp(uz,bz,lbeta);
    mpz_fdiv_q_2exp(bz,bz,lbeta);
    mpz_sub(tz,tz,uz);
    mpz_fdiv_r(tz,tz,qz);
    for(j=0;j<L;j++) {
      r->limbs[j].c[i] = mpz_get_ui(tz) & 0x3FFF;
      mpz_fdiv_q_2exp(tz,tz,14);
    }
  }
}

static void polz_kroneckermul_extension(polz *r, const polz *a, const polz *b, size_t deg) {
  ssize_t i,j,k;
  mpz_t az,bz,qz,tz,uz;
  const mp_bitcnt_t lbeta = 2*LOGQ+12;  // good for log(deg) <= 6

  mpz_inits(az,bz,qz,tz,uz,NULL);

  for(i=L-1;i>=0;i--) {
    mpz_mul_2exp(qz,qz,14);
    mpz_add_ui(qz,qz,modulus.q.limbs[i]);
  }

  for(i=N-1;i>=0;i--) {
    for(j=deg-1;j>=0;j--) {
      for(k=L-1;k>=0;k--) {
        mpz_add_ui(az,az,a[j].limbs[k].c[i]);
        mpz_add_ui(bz,bz,b[j].limbs[k].c[i]);
        if(k) {
          mpz_mul_2exp(az,az,14);
          mpz_mul_2exp(bz,bz,14);
        }
      }
      if(j) {
        mpz_mul_2exp(az,az,lbeta-14*(L-1));
        mpz_mul_2exp(bz,bz,lbeta-14*(L-1));
      }
    }
    if(i) {
      mpz_mul_2exp(az,az,lbeta-14*(L-1));
      mpz_mul_2exp(bz,bz,lbeta-14*(L-1));
    }
  }

  mpz_mul(az,az,bz);
  mpz_fdiv_q_2exp(bz,az,N*deg*lbeta);

  for(i=0;i<N;i++) {
    for(j=0;j<(ssize_t)deg;j++) {
      mpz_fdiv_r_2exp(tz,az,lbeta);
      mpz_fdiv_q_2exp(az,az,lbeta);
      mpz_fdiv_r_2exp(uz,bz,lbeta);
      mpz_fdiv_q_2exp(bz,bz,lbeta);
      mpz_sub(tz,tz,uz);
      mpz_fdiv_r(tz,tz,qz);
      for(k=0;k<L;k++) {
        r[j].limbs[k].c[i] = mpz_get_ui(tz) & 0x3FFF;
        mpz_fdiv_q_2exp(tz,tz,14);
      }
    }
  }
}

static void test_mul_extension(size_t deg) {
  polz a[deg], b[deg];
  polx ax[deg], bx[deg], cx[deg];
  poly t[deg];
  __attribute__((aligned(16)))
  uint8_t seed[16];
  uint64_t nonce = 0;

  randombytes(seed,16);
  polzvec_almostuniform(a,deg,seed,nonce++);
  polyvec_ternary(t,deg,seed,nonce++);
  polzvec_frompolyvec(b,t,deg);
  polzvec_reduce(b,deg);  // TODO: needed?
  polzvec_caddq(b,deg);

  polzvec_topolxvec(ax,a,deg);
  polxvec_frompolyvec(bx,t,deg);

  polz_kroneckermul_extension(a,a,b,deg);
  polxvec_mul_extension(cx,ax,bx,deg,deg,1);

  polzvec_frompolxvec(b,cx,deg);
  polzvec_sub(a,a,b,deg);
  polzvec_center(a,deg);
  if(!polzvec_iszero(a,deg))
    fprintf(stderr,"ERROR in polxvec_mul_extension()\n");

  polx alpha[deg], tmp0[1], tmp1[1];
  polxvec_quarternary(alpha,deg,seed,nonce++);
  polxvec_sprod(tmp0,alpha,cx,deg);
  polxvec_setzero(cx,deg);
  polxvec_collaps_add_extension(cx,alpha,ax,deg,deg,1);
  polxvec_sprod(tmp1,cx,bx,deg);
  polx_sub(tmp0,tmp0,tmp1);
  if(!polx_iszero(tmp0))
    fprintf(stderr,"ERROR in polxvec_collaps_add_extension()\n");
}

int main(void) {
  size_t i;
  size_t t,e;
  unsigned long long tsc[21], overhead;
  __attribute__((aligned(16)))
  uint8_t buf[N*QBYTES];
  uint64_t nonce = 0;
  polz a,b,c,d;
  polx x;
  poly y[LOGQ];
  double l;

  randombytes(buf,16);
  overhead = cpucycles_overhead();

  polzvec_uniform(&a,1,buf,nonce++);
  l = polzvec_norm(&a,1);
  printf("norm uniform: %.2g (rel dist to expected: %.2f)\n",l,fabs(1-l/ldexp(sqrt(N/12.0),LOGQ)));
  polzvec_almostuniform(&b,1,buf,nonce++);
  l = polzvec_norm(&b,1);
  printf("norm almost uniform: %.2g (rel dist to expected: %.2f)\n",l,fabs(1-l/ldexp(sqrt(N/12.0),LOGQ)));

  polz_bitpack(buf,&a);
  polz_bitunpack(&c,buf);
  polz_sub(&c,&c,&a);
  polz_reduce(&c);
  polz_center(&c);
  if(!polz_iszero(&c))
      fprintf(stderr,"ERROR in polz_bitpack/bitunpack\n");

  polz_topolx(&x,&a);
  polz_frompolx(&c,&x);
  polz_sub(&c,&c,&a);
  polz_reduce(&c);
  polz_center(&c);
  if(!polz_iszero(&c))
      fprintf(stderr,"ERROR in polz_to/frompolx\n");

  for(i=0;i<21;i++) {
    tsc[i] = cpucycles();
    polz_kroneckermul(&d,&a,&b);
  }
  for(i=0;i<20;i++)
    printf("polz_kroneckermul:  %2lu: %8lld\n", i, tsc[i+1] - tsc[i] - overhead);

  for(i=0;i<21;i++) {
    tsc[i] = cpucycles();
    polz_mul(&c,&a,&b);
  }
  for(i=0;i<20;i++)
    printf("polz:  %2lu: %8lld\n", i, tsc[i+1] - tsc[i] - overhead);

  polz_sub(&d,&d,&c);
  polz_reduce(&c);
  polz_center(&d);
  if(!polz_iszero(&d))
      fprintf(stderr,"ERROR in polz_mul\n");

  test_mul_extension(16);

  e = 7;
  t = (LOGQ+e-1)/e;
  polz_center(&a);
  polz_split(y,&c,&a,e);
  polz_slli(&c,&c,e);
  polz_frompoly(&d,y);
  polz_add(&c,&c,&d);
  polz_sub(&c,&c,&a);
  polz_reduce(&c);
  polz_center(&c);
  if(!polz_iszero(&c))
      fprintf(stderr,"ERROR in polz_split\n");

  polz_decompose(y,&a,1,t,e);
  polz_frompoly(&c,&y[t-1]);
  for(i=1;i<t;i++) {
    polz_slli(&c,&c,e);
    polz_frompoly(&d,&y[t-1-i]);
    polz_add(&c,&c,&d);
  }
  polz_sub(&c,&c,&a);
  polz_reduce(&c);
  polz_center(&c);
  if(!polz_iszero(&c))
      fprintf(stderr,"ERROR in polz_decompose / polz_slli\n");

  polz_decompose(y,&a,1,t,e);
  polx_frompoly(&x,&y[0]);
  for(i=1;i<t;i++) {
    polx tmp;
    polx_frompoly(&tmp,&y[i]);
    polx_scale_add(&x,&tmp,(int64_t)1 << e*i);
  }
  polz_frompolx(&c,&x);
  polz_sub(&c,&c,&a);
  polz_reduce(&c);
  polz_center(&c);
  if(!polz_iszero(&c))
      fprintf(stderr,"ERROR in polz_decompose\n");

  polz_decompose(y,&a,1,t,e);
  polz_reconstruct(&c,y,1,t,e);
  polz_sub(&c,&c,&a);
  polz_reduce(&c);
  polz_center(&c);
  if(!polz_iszero(&c))
      fprintf(stderr,"ERROR in polz_reconstruct\n");

  return 0;
}

