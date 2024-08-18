#include <stdio.h>
#include <stdlib.h>
#include "randombytes.h"
#include "fips202.h"
#include "dachshund.h"
#include "pack.h"

static void prepare_linear(smplstmnt *st, witness *wt) {
  size_t i,j,l;
  __attribute__((aligned(16)))
  uint8_t seed[16];
  uint64_t nonce = 0;
  shake128incctx shakectx;
  sparsecnst *cnst;
  polx *buf;

  size_t r = 2;
  size_t n[r];
  for(i=0;i<r;i++)
    n[i] = 1<<16;
  size_t k = r;
  size_t deg = 8;
  size_t deg2 = 8; // next2power
  size_t betasq[r];
  for(i=0;i<r;i++)
    betasq[i] = 1.15*10/16*n[i]*N;

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
  init_smplstmnt_raw(st,r,n,betasq,k);
  for(i=0;i<k;i++) {
    cnst = &st->cnst[i];
    l = extlen(n[i],deg);
    buf = init_sparsecnst_half(cnst,r,1,l,deg,0,0);

    cnst->idx[0] = i;
    cnst->off[0] = 0;
    cnst->len[0] = n[i];
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
  commitment com = {};
  smplstmnt st0 = {};
  statement st1 = {}, st2 = {};
  proof pi0 = {}, pi1 = {};
  witness wt0 = {}, wt1 = {}, wt2 = {};
  double size = 0;

  printf("Testing Dachshund followed by one Labrador\n\n");

  prepare_linear(&st0,&wt0);
  print_smplstmnt_pp(&st0);
  ret = simple_verify(&st0,&wt0);
  if(ret) {
    fprintf(stderr,"ERROR: verification of prepare_linear failed: %d\n",ret);
    goto end;
  }

  simple_prove(&st1,&wt1,&pi0,&com,&st0,&wt0,0);
  free_witness(&wt0);
  size += print_proof_pp(&pi0);
  print_statement_pp(&st1);
  ret = verify(&st1,&wt1);
  if(ret) {
    fprintf(stderr,"ERROR: verification of Dachshund failed: %d\n",ret);
    goto end;
  }

  free_statement(&st1);
  ret = simple_reduce(&st1,&pi0,&com,&st0);
  free_smplstmnt(&st0);
  if(ret) {
    fprintf(stderr,"ERROR: simple_reduce failed: %d\n",ret);
    goto end;
  }
  ret = verify(&st1,&wt1);
  if(ret) {
    fprintf(stderr,"ERROR: verification of simple_reduce failed: %d\n",ret);
    goto end;
  }

  prove(&st2,&wt2,&pi1,&st1,&wt1,0);
  free_witness(&wt1);
  size += print_proof_pp(&pi1);
  print_statement_pp(&st2);
  ret = verify(&st2,&wt2);
  if(ret) {
    fprintf(stderr,"ERROR: verification of Labrador failed: %d\n",ret);
    goto end;
  }

  free_statement(&st2);
  ret = reduce(&st2,&pi1,&st1);
  free_statement(&st1);
  if(ret) {
    fprintf(stderr,"ERROR: reduce failed: %d\n",ret);
    goto end;
  }
  ret = verify(&st2,&wt2);
  if(ret) {
    fprintf(stderr,"ERROR: verification of reduce failed: %d\n",ret);
    goto end;
  }

  size += print_witness_pp(&wt2);
  printf("Total proof size: %.2f KB\n",size);
  printf("\n");

end:
  free_commitment(&com);
  free_smplstmnt(&st0);
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
  smplstmnt st = {};
  witness wt = {};
  commitment com = {};
  composite p = {};

  printf("Testing Dachshund Composite\n\n");

  prepare_linear(&st,&wt);
  print_smplstmnt_pp(&st);
  ret = simple_verify(&st,&wt);
  if(ret) {
    fprintf(stderr,"ERROR: verification of prepare_linear failed: %d\n",ret);
    goto end;
  }

  composite_prove_simple(&p,&com,&st,&wt);;
  ret = composite_verify_simple(&p,&com,&st);
  if(ret)
    fprintf(stderr,"ERROR: verify_composite_simple failed: %d\n",ret);

end:
  free_smplstmnt(&st);
  free_commitment(&com);
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
