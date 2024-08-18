#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "malloc.h"
#include "randombytes.h"
#include "labrador.h"
#include "chihuahua.h"
#include "pack.h"
#include "greyhound.h"

static int test_polcom(size_t len) {
  int ret;
  int64_t x,y;
  polz *s;
  polcomctx ctx;
  polcomprf pi;
  prncplstmnt st;
  witness wt;
  __attribute__((aligned(16)))
  uint8_t seed[16];

  printf("Testing Greyhound polynomial commitment scheme\n\n");
  randombytes(seed,16);
  s = _aligned_alloc(64,len*sizeof(polz));
  polzvec_almostuniform(s,len,seed,0);
  polzvec_center(s,len);

  x = 43;
  y = polzvec_eval(s,len,x);

  polcom_commit(&ctx,s,len);
  print_polcomctx_pp(&ctx);

  polcom_eval(&wt,&pi,&ctx,x,y);
  free(s);
  free_polcomctx(&ctx);
  print_polcomprf_pp(&pi);

  ret = polcom_reduce(&st,&pi);
  if(ret) {
    printf("ERROR: Reduction to Chihuahua statement failed: %d\n",ret);
    goto end;
  }
  free_polcomprf(&pi);
  print_prncplstmnt_pp(&st);
  print_witness_pp(&wt);

  ret = principle_verify(&st,&wt);
  if(ret) {
    printf("ERROR: Verification of Chihuahua statement failed: %d\n",ret);
    goto end;
  }

end:
  free_polcomprf(&pi);
  free_prncplstmnt(&st);
  free_witness(&wt);
  return ret;
}

static int test_pack(size_t len) {
  int ret;
  int64_t x,y;
  clock_t t;
  polz *s;
  polcomctx ctx = {};
  polcomprf pi = {};
  composite p = {};
  __attribute__((aligned(16)))
  uint8_t seed[16];

  printf("Testing Greyhound Pack for degree 2^%.2g\n\n",log2(len << 6));
  randombytes(seed,16);
  s = _aligned_alloc(64,len*sizeof(polz));
  polzvec_almostuniform(s,len,seed,0);
  polzvec_center(s,len);

  x = 43;
  y = polzvec_eval(s,len,x);

  t = clock();
  polcom_commit(&ctx,s,len);
  t = clock() - t;
  printf("Greyhound Pack commit time: %.4fs\n\n",(double)t/CLOCKS_PER_SEC);
  print_polcomctx_pp(&ctx);

  composite_prove_polcom(&p,&pi,&ctx,x,y);
  ret = composite_verify_polcom(&p,&pi);
  if(ret)
    printf("ERROR: verify_composite_polcom failed: %d\n",ret);

  free(s);
  free_polcomctx(&ctx);
  free_polcomprf(&pi);
  free_composite(&p);
  return ret;
}

int main(void) {
  int ret;
  size_t len;

  len = 1 << 20;

  ret = test_polcom(len);
  if(ret) goto end;
  ret = test_pack(len);
  if(ret) goto end;

end:
  free_comkey();
  return ret;
}
