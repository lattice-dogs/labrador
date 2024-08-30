#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "randombytes.h"
#include "fips202.h"
#include "chihuahua.h"
#include "pack.h"

static void prepare_linear(prncplstmnt *st, witness *wt) {
  size_t i,j,l;
  __attribute__((aligned(16)))
  uint8_t seed[16];
  uint64_t nonce = 0;
  shake128incctx shakectx;
  sparsecnst *cnst;
  polx *buf;

  size_t r = 1;
  size_t n[r];
  for(i=0;i<r;i++)
    n[i] = 1<<11;
  size_t k = 2;
  size_t deg = 8;
  size_t deg2 = 8; // next2power
  size_t betasq = 0;
  for(i=0;i<r;i++)
    betasq += 1.15*10/16*n[i]*N;

  __attribute__((aligned(16)))
  uint8_t hashbuf[deg2*N*QBYTES];
  polx *sx[r];
  polz t[deg2];

  randombytes(seed,sizeof(seed));
  init_witness_raw(wt,r,n);
  for(i=0;i<r;i++) {
    polyvec_ternary(wt->s[i],wt->n[i],seed,nonce++);
    wt->normsq[i] = polyvec_sprodz(wt->s[i],wt->s[i],wt->n[i]);
  }

  *sx = NULL;
  shake128_inc_init(&shakectx);
  init_prncplstmnt_raw(st,r,n,betasq,k,0);
  for(i=0;i<k;i++) {
    cnst = &st->cnst[i];
    l = extlen(n[0],deg);
    buf = init_sparsecnst_half(cnst,r,1,l,deg,0,0);

    cnst->idx[0] = 0;
    cnst->off[0] = 0;
    cnst->len[0] = n[0];
    cnst->mult[0] = 1;
    cnst->phi[0] = buf;

    for(j=0;j<l;j+=deg2) {
      polzvec_almostuniform(t,deg2,seed,nonce++);
      polzvec_bitpack(hashbuf,t,deg2);
      shake128_inc_absorb(&shakectx,hashbuf,deg2*N*QBYTES);
      polzvec_topolxvec(&cnst->phi[0][j],t,deg2);
    }

    sparsecnst_eval(cnst->b,cnst,sx,wt);
    polzvec_frompolxvec(t,cnst->b,deg);
    polzvec_bitpack(hashbuf,t,deg);
    shake128_inc_absorb(&shakectx,hashbuf,deg*N*QBYTES);
  }

  free(*sx);
  shake128_inc_finalize(&shakectx);
  shake128_inc_squeeze(st->h,16,&shakectx);
}

static int test_twolayer() {
  int ret;
  prncplstmnt st0 = {};
  statement st1 = {}, st2 = {};
  proof pi0 = {}, pi1 = {};
  witness wt0 = {}, wt1 = {}, wt2 = {};
  double size = 0;

  printf("Testing Chihuahua followed by one Labrador\n\n");

  prepare_linear(&st0,&wt0);
  print_prncplstmnt_pp(&st0);
  ret = principle_verify(&st0,&wt0);
  if(ret) {
    fprintf(stderr,"ERROR: Verification of prepare_linear failed: %d\n",ret);
    goto end;
  }

  ret = principle_prove(&st1,&wt1,&pi0,&st0,&wt0,0);
  if(ret) {
    fprintf(stderr,"ERROR: Chihuahua proof failed: %d\n",ret);
    goto end;
  }
  free_witness(&wt0);
  size += print_proof_pp(&pi0);
  print_statement_pp(&st1);
  ret = verify(&st1,&wt1);
  if(ret) {
    fprintf(stderr,"ERROR: Chihuahua verification failed: %d\n",ret);
    goto end;
  }

  free_statement(&st1);
  ret = principle_reduce(&st1,&pi0,&st0);
  free_prncplstmnt(&st0);
  if(ret) {
    fprintf(stderr,"ERROR: Chihuahua reduction failed: %d\n",ret);
    goto end;
  }
  ret = verify(&st1,&wt1);
  if(ret) {
    fprintf(stderr,"ERROR: Verification of chihuahua reduction failed: %d\n",ret);
    goto end;
  }

  ret = prove(&st2,&wt2,&pi1,&st1,&wt1,0);
  if(ret) {
    fprintf(stderr,"ERROR: Labrador proof failed: %d\n",ret);
    goto end;
  }
  free_witness(&wt1);
  size += print_proof_pp(&pi1);
  print_statement_pp(&st2);
  ret = verify(&st2,&wt2);
  if(ret) {
    fprintf(stderr,"ERROR: Labrador verification failed: %d\n",ret);
    goto end;
  }

  free_statement(&st2);
  ret = reduce(&st2,&pi1,&st1);
  free_statement(&st1);
  if(ret) {
    fprintf(stderr,"ERROR: Labrador reduction failed: %d\n",ret);
    goto end;
  }
  ret = verify(&st2,&wt2);
  if(ret) {
    fprintf(stderr,"ERROR: Verification of Labrador reduction failed: %d\n",ret);
    goto end;
  }

  size += print_witness_pp(&wt2);
  printf("Total proof size: %.2f KB\n",size);
  printf("\n");

end:
  free_prncplstmnt(&st0);
  free_statement(&st1);
  free_statement(&st2);
  free_proof(&pi0);
  free_proof(&pi1);
  free_witness(&wt0);
  free_witness(&wt1);
  free_witness(&wt2);
  return ret;
}

static int test_pack() {
  int ret;
  prncplstmnt st = {};
  witness wt = {};
  composite p = {};

  printf("Testing Chihuahua Composite\n\n");

  prepare_linear(&st,&wt);
  print_prncplstmnt_pp(&st);
  ret = principle_verify(&st,&wt);
  if(ret) {
    fprintf(stderr,"ERROR: Verification of prepare_linear failed: %d\n",ret);
    goto end;
  }

  ret = composite_prove_principle(&p,&st,&wt);
  if(ret) {
    fprintf(stderr,"ERROR: Chihuahua composite proof failed: %d\n",ret);
    goto end;
  }
  ret = composite_verify_principle(&p,&st);
  if(ret) {
    fprintf(stderr,"ERROR: Chihuahua composite verifaction failed: %d\n",ret);
    goto end;
  }

end:
  free_prncplstmnt(&st);
  free_composite(&p);
  free_witness(&wt);
  return ret;
}

int main(void) {
  int ret;

  ret = test_twolayer();
  if(ret) goto end;
  ret = test_pack();
  if(ret) goto end;

end:
  free_comkey();
  return ret;
}
