#include <stdlib.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include "malloc.h"
#include "labrador.h"
#include "chihuahua.h"
#include "dachshund.h"
#include "greyhound.h"
#include "pack.h"

void free_composite(composite *p) {
  size_t i;

  free_witness(&p->owt);
  for(i=0;i<p->l;i++) {
    free_proof(p->pi[i]);
    free(p->pi[i]);
    p->pi[i] = NULL;
  }
  p->l = 0;
}

static int composite_prove(composite *p, statement *tst, witness *twt, double *twtsize) {
  int ret;
  size_t i = 0;
  double pisize;

  while(p->l < 16) {
    p->pi[p->l] = _malloc(sizeof(proof));
    ret = prove(&tst[i^1],&twt[i^1],p->pi[p->l],&tst[i],&twt[i],0);
    if(ret) return ret;
    pisize = print_proof_pp(p->pi[p->l]);
    print_statement_pp(&tst[i^1]);
    twtsize[i^1] = print_witness_pp(&twt[i^1]);
    if(pisize + twtsize[i^1] >= twtsize[i]) {
      free_proof(p->pi[p->l]);
      free_statement(&tst[i^1]);
      free_witness(&twt[i^1]);
      break;
    }

    free_statement(&tst[i]);
    free_witness(&twt[i]);
    p->size += pisize;
    p->l += 1;
    i ^= 1;
  }

  if(p->l < 16) {
    ret = prove(&tst[i^1],&twt[i^1],p->pi[p->l],&tst[i],&twt[i],1);
    if(ret) return ret;
    pisize = print_proof_pp(p->pi[p->l]);
    print_statement_pp(&tst[i^1]);
    twtsize[i^1] = print_witness_pp(&twt[i^1]);
    if(pisize + twtsize[i^1] >= twtsize[i]) {
      free_proof(p->pi[p->l]);
      free(p->pi[p->l]);
    }
    else {
      p->size += pisize;
      p->l += 1;
      i ^= 1;
    }

    free_statement(&tst[i^1]);
    free_witness(&twt[i^1]);
  }

  free_statement(&tst[i]);
  p->owt = twt[i];
  p->size += twtsize[i];
  return 0;
}

int composite_prove_principle(composite *p, const prncplstmnt *st, const witness *wt) {
  int ret;
  statement tst[2] = {};
  witness twt[2] = {};
  double twtsize[2];
  clock_t t;

  p->l = 0;
  memset(&p->owt,0,sizeof(witness));
  p->pi[p->l] = _malloc(sizeof(proof));

  t = clock();
  ret = principle_prove(tst,twt,p->pi[p->l],st,wt,0);
  if(ret)
    goto err;
  p->size = print_proof_pp(p->pi[p->l]);
  print_statement_pp(tst);
  twtsize[0] = print_witness_pp(twt);
  p->l += 1;
  ret = composite_prove(p,tst,twt,twtsize);
  if(ret) {
    ret += 10;
    goto err;
  }
  t = clock() - t;

  printf("Commitment key length: %zu\n",comkey_len);
  printf("Chihuahua Pack members: %zu\n",p->l);
  printf("Chihuahua Pack size: %.2f KB\n",p->size);
  printf("Chihuahua Pack proving time: %.4fs\n",(double)t/CLOCKS_PER_SEC);
  return 0;

err:
  free_statement(&tst[0]);
  free_statement(&tst[1]);
  free_witness(&twt[0]);
  free_witness(&twt[1]);
  free(p->pi[p->l]);
  p->pi[p->l] = NULL;
  free_composite(p);
  return ret;
}

int composite_prove_simple(composite *p, commitment *com, const smplstmnt *st, const witness *wt) {
  int ret;
  statement tst[2] = {};
  witness twt[2] = {};
  double twtsize[2];
  clock_t t;

  p->l = 0;
  memset(&p->owt,0,sizeof(witness));
  p->pi[p->l] = _malloc(sizeof(proof));

  t = clock();
  ret = simple_prove(tst,twt,p->pi[p->l],com,st,wt,0);
  if(ret)
    goto err;
  p->size = print_proof_pp(p->pi[p->l]);
  print_statement_pp(tst);
  twtsize[0] = print_witness_pp(twt);
  p->l += 1;

  ret = composite_prove(p,tst,twt,twtsize);
  if(ret) {
    ret += 10;
    goto err;
  }
  t = clock() - t;

  printf("Commitment key length: %zu\n",comkey_len);
  printf("Dachshund Pack members: %zu\n",p->l);
  printf("Dachshund Pack size: %.2f KB\n",p->size);
  printf("Dachshund Pack proving time: %.4fs\n",(double)t/CLOCKS_PER_SEC);
  return 0;

err:
  free_statement(&tst[0]);
  free_statement(&tst[1]);
  free_witness(&twt[0]);
  free_witness(&twt[1]);
  free(p->pi[p->l]);
  p->pi[p->l] = NULL;
  free_composite(p);
  return ret;
}

int composite_prove_polcom(composite *p, polcomprf *ppi, polcomctx *ctx, uint32_t x, uint32_t y) {
  int ret;
  prncplstmnt tst0[1] = {};
  statement tst[2] = {};
  witness twt[2] = {};
  double twtsize[2];
  clock_t t;

  p->l = 0;
  memset(&p->owt,0,sizeof(witness));
  p->pi[p->l] = _malloc(sizeof(proof));

  t = clock();
  polcom_eval(&twt[1],ppi,ctx,x,y);
  ret = polcom_reduce(tst0,ppi);
  if(ret)
    goto err;
  p->size = print_polcomprf_pp(ppi);
  print_prncplstmnt_pp(tst0);
  twtsize[1] = print_witness_pp(&twt[1]);

  ret = principle_prove(tst,twt,p->pi[p->l],tst0,&twt[1],0);
  if(ret) {
    ret += 10;
    goto err;
  }
  p->size += print_proof_pp(p->pi[p->l]);
  print_statement_pp(tst);
  twtsize[0] = print_witness_pp(twt);
  free_prncplstmnt(tst0);
  free_witness(&twt[1]);
  p->l += 1;

  ret = composite_prove(p,tst,twt,twtsize);
  if(ret) {
    ret += 20;
    goto err;
  }
  t = clock() - t;

  printf("Commitment key length: %zu\n",comkey_len);
  printf("Greyhound Pack members: %zu\n",p->l);
  printf("Greyhound Pack size: %.2f KB\n",p->size);
  printf("Greyhound Pack proving time: %.4fs\n",(double)t/CLOCKS_PER_SEC);
  return 0;

err:
  free_prncplstmnt(tst0);
  free_statement(&tst[0]);
  free_statement(&tst[1]);
  free_witness(&twt[0]);
  free_witness(&twt[1]);
  free(p->pi[p->l]);
  p->pi[p->l] = NULL;
  free_composite(p);
  return ret;
}

static int composite_verify(const composite *p, statement *tst) {
  size_t i,j;
  int ret;

  i = 0;
  for(j=1;j<p->l;j++) {
    ret = reduce(&tst[i^1],p->pi[j],&tst[i]);
    free_statement(&tst[i]);
    i ^= 1;
    if(ret)  // projection too long or commitments not secure (1/2/3)
      return ret + 4*j;
  }

  ret = verify(&tst[i],&p->owt);
  if(ret)
    return ret + 100;

  return 0;
}

int composite_verify_principle(const composite *p, const prncplstmnt *st) {
  int ret;
  statement tst[2] = {};
  clock_t t;

  t = clock();
  ret = principle_reduce(tst,p->pi[0],st);
  if(ret)  // projection too long or commitments not secure (1/2/3)
    goto err;

  ret = composite_verify(p,tst);
  if(ret) {
    ret += 10;
    goto err;
  }
  t = clock() - t;
  printf("Chihuahua Pack verification time: %.4fs\n",(double)t/CLOCKS_PER_SEC);
  return 0;

err:
  free_statement(&tst[0]);
  free_statement(&tst[1]);
  return ret;
}

int composite_verify_simple(const composite *p, const commitment *com, const smplstmnt *st) {
  int ret;
  statement tst[2] = {};
  clock_t t;

  t = clock();
  ret = simple_reduce(tst,p->pi[0],com,st);
  if(ret)
    goto err;

  ret = composite_verify(p,tst);
  if(ret) {
    ret += 10;
    goto err;
  }
  t = clock() - t;
  printf("Dachshund Pack verification time: %.4fs\n",(double)t/CLOCKS_PER_SEC);
  return 0;

err:
  free_statement(&tst[0]);
  free_statement(&tst[1]);
  return ret;
}

int composite_verify_polcom(const composite *p, const polcomprf *ppi) {
  int ret;
  prncplstmnt tst0[1] = {};
  statement tst[2] = {};
  clock_t t;

  t = clock();
  ret = polcom_reduce(tst0,ppi);
  if(ret)
    goto err;

  ret = principle_reduce(tst,p->pi[0],tst0);
  free_prncplstmnt(tst0);
  if(ret) {
    ret += 10;
    goto err;
  }

  ret = composite_verify(p,tst);
  if(ret) {
    ret += 20;
    goto err;
  }
  t = clock() - t;
  printf("Greyhound Pack verification time: %.4fs\n",(double)t/CLOCKS_PER_SEC);
  return 0;

err:
  free_prncplstmnt(tst0);
  free_statement(&tst[0]);
  free_statement(&tst[1]);
  return ret;
}
