#include <stdlib.h>
#include <stddef.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <stdio.h>
#include "malloc.h"
#include "fips202.h"
#include "aesctr.h"
#include "data.h"
#include "poly.h"
#include "polx.h"
#include "polz.h"
#include "sparsemat.h"
#include "jlproj.h"
#include "labrador.h"
#include "chihuahua.h"

polx *init_sparsecnst_half(sparsecnst *cnst, size_t r, size_t nz, size_t buflen, size_t deg,
                           int quadratic, int homogeneous)
{
  void *buf;

  cnst->deg = deg;
  cnst->nz = nz;
  cnst->a->len = 0;
  if(quadratic) {
    buf = _malloc((r*r+r)*sizeof(size_t));
    cnst->a->rows = (size_t*)buf;
    cnst->a->cols = (size_t*)buf + (r*r+r)/2;

    buf = _aligned_alloc(64,r*sizeof(polx));
    cnst->a->coeffs = (polx*)buf;
  }
  else {
    cnst->a->rows = NULL;
    cnst->a->cols = NULL;
    cnst->a->coeffs = NULL;
  }

  if(!homogeneous)
    buflen += MAX(1,deg);

  buf = _malloc(4*nz*sizeof(size_t) + nz*sizeof(polx*) + 63 + buflen*sizeof(polx));
  cnst->idx = (size_t*)buf;
  cnst->off = &cnst->idx[nz];
  cnst->len = &cnst->off[nz];
  cnst->mult = &cnst->len[nz];
  cnst->phi = (polx**)&cnst->mult[nz];
  buf = (void*)(((uintptr_t)&cnst->phi[nz] + 63) & -64ULL);

  if(homogeneous)
    cnst->b = NULL;
  else {
    cnst->b = buf;
    buf = (polx*)buf + MAX(1,deg);
  }

  return buf;
}

void init_sparsecnst_raw(sparsecnst *cnst, size_t r, size_t nz, const size_t idx[nz], const size_t n[nz], size_t deg,
                         int quadratic, int homogeneous)
{
  size_t i;
  size_t buflen;
  size_t elen[nz];
  polx *buf;

  buflen = 0;
  for(i=0;i<nz;i++) {
    elen[i] = extlen(n[i]*deg,deg);
    buflen += elen[i];
  }

  buf = init_sparsecnst_half(cnst,r,nz,buflen,deg,quadratic,homogeneous);
  for(i=0;i<nz;i++) {
    cnst->idx[i] = idx[i];
    cnst->off[i] = 0;
    cnst->len[i] = n[i]*deg;
    cnst->mult[i] = 1;
    cnst->phi[i] = (polx*)buf;
    buf += elen[i];
  }
}

int set_sparsecnst_raw(sparsecnst *cnst, uint8_t h[16], size_t nz, const size_t idx[nz],
                       const size_t n[nz], size_t deg, int64_t *phi, int64_t *b)
{
  size_t j,k;
  polz t[deg];
  __attribute__((aligned(16)))
  uint8_t hashbuf[deg*N*QBYTES];
  shake128incctx shakectx;

  if(nz != cnst->nz)
    return 1;
  if(deg != cnst->deg)
    return 2;
  if((b && !cnst->b) || (!b && cnst->b))
    return 3;
  for(j=0;j<nz;j++) {
    if(idx[j] != cnst->idx[j])
      return 4;
    if(cnst->off[j])
      return 5;
    if(cnst->mult[j] != 1)
      return 6;
  }

  cnst->a->len = 0;

  shake128_inc_init(&shakectx);
  shake128_inc_absorb(&shakectx,h,16);

  if(b) {
    polzvec_fromint64vec(t,1,deg,b);
    polzvec_topolxvec(cnst->b,t,deg);
    polzvec_bitpack(hashbuf,t,deg);
    shake128_inc_absorb(&shakectx,hashbuf,deg*N*QBYTES);
  }

  for(j=0;j<nz;j++) {
    for(k=0;k<n[j];k++) {
      polzvec_fromint64vec(t,1,deg,phi);
      polzvec_topolxvec(&cnst->phi[j][k*deg],t,deg);
      polzvec_bitpack(hashbuf,t,deg);
      shake128_inc_absorb(&shakectx,hashbuf,deg*N*QBYTES);
      phi += deg*N;
    }
  }

  shake128_inc_finalize(&shakectx);
  shake128_inc_squeeze(h,16,&shakectx);
  return 0;
}

void free_sparsecnst(sparsecnst *cnst) {
  free(cnst->idx);
  free(cnst->a->rows);
  free(cnst->a->coeffs);
  cnst->idx = cnst->a->rows = NULL;
  cnst->a->coeffs = NULL;
}

void sparsecnst_eval(polx *b, const sparsecnst *cnst, polx *sx[], const witness *wt) {
  const size_t r = wt->r;
  const size_t *n = wt->n;

  size_t i,j,k;
  size_t deg2 = MAX(1,cnst->deg);
  polx t[deg2];

  if(!*sx) {
    k = 0;
    for(i=0;i<r;i++)
      k += n[i];
    sx[0] = _aligned_alloc(64,k*sizeof(polx));
    for(i=1;i<r;i++)
      sx[i] = &sx[i-1][n[i-1]];
    for(i=0;i<r;i++)
      polxvec_frompolyvec(sx[i],wt->s[i],n[i]);
  }

  polxvec_setzero(b,deg2);

  /* quadratic term */
  for(i=0;i<cnst->a->len;i++) {
    j = cnst->a->rows[i];
    k = cnst->a->cols[i];
    polxvec_sprod(t,sx[j],sx[k],MIN(n[j],n[k]));  // TODO: store?
    polx_mul(t,t,&cnst->a->coeffs[i]);
    polx_add(b,b,t);
  }

  /* linear term */
  for(i=0;i<cnst->nz;i++) {
    j = cnst->idx[i];
    k = cnst->off[i];
    polxvec_mul_extension(t,cnst->phi[i],&sx[j][k],cnst->len[i],deg2/cnst->mult[i],cnst->mult[i]);
    polxvec_add(b,b,t,deg2);
  }
}

int sparsecnst_check(const sparsecnst *cnst, polx *sx[], const witness *wt) {
  int ret;
  size_t deg2 = MAX(1,cnst->deg);
  polx b[deg2];

  sparsecnst_eval(b,cnst,sx,wt);
  if(cnst->b)
    polxvec_sub(b,b,cnst->b,deg2);

  if(cnst->deg)
    ret = polxvec_iszero(b,deg2);
  else
    ret = polx_iszero_constcoeff(b);

  return ret;
}

void init_prncplstmnt_raw(prncplstmnt *st, size_t r, const size_t n[r],
                          uint64_t betasq, size_t k, int quadratic)
{
  size_t i;
  void *buf;

  buf = _malloc(r*sizeof(size_t) + k*sizeof(sparsecnst));
  st->n = (size_t*)buf;
  st->cnst = (sparsecnst*)&st->n[r];

  st->r = r;
  st->k = k;
  st->quadratic = quadratic;
  st->betasq = betasq;
  for(i=0;i<r;i++)
    st->n[i] = n[i];
  for(i=0;i<k;i++)
    st->cnst[i].idx = NULL;

  memset(st->h,0,16);
}

int set_prncplstmnt_lincnst_raw(prncplstmnt *st, size_t i, size_t nz, const size_t idx[nz],
                                const size_t n[nz], size_t deg, int64_t *phi, int64_t *b)
{
  size_t j,k;

  if(i >= st->k)
    return 1;

  sparsecnst *cnst = &st->cnst[i];
  if(cnst->idx)
    return 2;
  for(j=0;j<nz;j++) {
    k = idx[j];
    if(k >= st->r)
      return 3;
    if(n[j]*deg != st->n[k])
      return 4;
  }

  init_sparsecnst_raw(cnst,st->r,nz,idx,n,deg,0,b == NULL);
  set_sparsecnst_raw(cnst,st->h,nz,idx,n,deg,phi,b);
  return 0;
}

void free_prncplstmnt(prncplstmnt *st) {
  size_t i;

  if(!st->n) return;
  for(i=0;i<st->k;i++)
    free_sparsecnst(&st->cnst[i]);
  free(st->n);  // one buffer for n,cnst
  st->n = NULL;
}

void print_prncplstmnt_pp(const prncplstmnt *st) {
  size_t i;

  printf("Chihuahua statement ");
  for(i=0;i<5;i++)
    printf("%02hhX",st->h[i]);
  printf(":\n");

  printf("  Witness multiplicity: %zu\n",st->r);
  printf("  Witness ranks: ");
  for(i=0;i<st->r;i++) {
    printf("%zu",st->n[i]);
    if(i<st->r-1) printf(", ");
  }
  printf("\n");

  printf("  Number of dot-product constraints: %zu\n",st->k);
  printf("  Norm constraint: %.2f\n",sqrt(st->betasq));
  printf("\n");
}

void collaps_sparsecnst(constraint *ocnst, statement *ost, const proof *pi, const sparsecnst *icnst, size_t k) {
  size_t i,j,u,v;
  const size_t n = ost->n;
  int64_t s;
  polx (*phi)[n];

  j = 0;
  for(i=0;i<k;i++)
    if(icnst[i].deg == 0)
      j += 1;

  __attribute__((aligned(16)))
  uint8_t hashbuf[16+j*QBYTES];

  shake128(hashbuf,sizeof(hashbuf),ost->h,16);
  memcpy(ost->h,hashbuf,16);

  for(i=0;i<k;i++) {
    if(icnst->deg) continue;

    s = 0;
    for(j=0;j<QBYTES;j++)
      s |= (int64_t)hashbuf[16+i*QBYTES+j] << 8*j;
    s &= (1LL << LOGQ) - 1;

    if(icnst->b)
      polx_scale_add(ocnst->b,icnst->b,s);

    /* Assumes sparse constraint is ordered */
    phi = (polx(*)[n])ocnst->phi;
    u = v = 0;
    for(j=0;j<icnst->nz;j++) {
      while(u < icnst->idx[j]) {
        phi += pi->nu[u];
        v = (pi->nu[u]) ? 0 : v + pi->n[u];
        u += 1;
      }
      polxvec_scale_add(&(*phi)[v+icnst->off[j]],icnst->phi[j],icnst->len[j],s);
    }

    sparsemat_scale_add(ocnst->a,icnst->a,s);
    icnst += 1;
  }
}

void aggregate_sparsecnst(statement *ost, const proof *pi, const sparsecnst *cnst, size_t k) {
  size_t i,j,u,v;
  __attribute__((aligned(16)))
  uint8_t hashbuf[32];
  polx (*phi)[ost->n];

  shake128(hashbuf,32,ost->h,16);
  memcpy(ost->h,hashbuf,16);

  for(i=0;i<k;i++) {
    if(cnst->deg == 0) continue;
    polx alpha[cnst->deg];
    polxvec_quarternary(alpha,cnst->deg,&hashbuf[16],i);

    if(cnst->b)
      polxvec_sprod_add(ost->cnst->b,alpha,cnst->b,cnst->deg);

    phi = (polx(*)[ost->n])ost->cnst->phi;
    u = v = 0;
    for(j=0;j<cnst->nz;j++) {
      while(u < cnst->idx[j]) {
        phi += pi->nu[u];
        v = (pi->nu[u]) ? 0 : v + pi->n[u];
        u += 1;
      }
      polxvec_collaps_add_extension(&(*phi)[v+cnst->off[j]],alpha,cnst->phi[j],cnst->len[j],
                                    cnst->deg/cnst->mult[j],cnst->mult[j]);
    }

    sparsemat_polx_mul_add(ost->cnst->a,alpha,cnst->a);
    cnst += 1;
  }

  ost->cnst->a->coeffs = realloc(ost->cnst->a->coeffs,ost->cnst->a->len*sizeof(polx));
}

void principle_prove(statement *ost, witness *owt, proof *pi, const prncplstmnt *ist, const witness *iwt, int tail) {
  size_t i;
  constraint cnst[1];

  init_proof(pi,iwt,ist->quadratic,tail);
  init_statement(ost,pi,ist->h);
  init_witness(owt,ost);
  printf("Predicted witness norm: %.2f\n\n",sqrt(pi->normsq));

  polx (*sx)[ost->n] = _aligned_alloc(64,ost->r*ost->n*sizeof(polx));
  commit(ost,owt,pi,sx,iwt);

  uint8_t (*jlmat)[ost->n][256*N/8] = _aligned_alloc(64,ost->r*ost->n*256*N/8);
  project(ost,pi,jlmat,iwt);

  init_constraint(cnst,ost);
  for(i=0;i<LIFTS;i++) {
    collaps_jlproj(cnst,ost,pi,jlmat);
    collaps_sparsecnst(cnst,ost,pi,ist->cnst,ist->k);
    lift_aggregate_zqcnst(ost,pi,i,cnst,sx);
  }
  free(jlmat);
  free_constraint(cnst);

  aggregate_sparsecnst(ost,pi,ist->cnst,ist->k);
  amortize(ost,owt,pi,sx);
  free(sx);

  polx_refresh(ost->cnst->b);
  polxvec_refresh(ost->cnst->phi,ost->n);
}

int principle_reduce(statement *ost, const proof *pi, const prncplstmnt *ist) {
  size_t i;
  int ret;
  uint8_t (*jlmat)[ost->n][256*N/8];
  constraint cnst[1];

  init_statement(ost,pi,ist->h);
  jlmat = _aligned_alloc(64,ost->r*ost->n*256*N/8);

  reduce_commit(ost,pi);
  ret = reduce_project(ost,jlmat,pi,pi->r,ist->betasq);
  if(ret) goto err;

  init_constraint(cnst,ost);
  for(i=0;i<LIFTS;i++) {
    collaps_jlproj(cnst,ost,pi,jlmat);
    collaps_sparsecnst(cnst,ost,pi,ist->cnst,ist->k);
    ret = reduce_lift_aggregate_zqcnst(ost,pi,i,cnst);
    if(ret) {
      ret = 2+i;
      goto err;
    }
  }
  free_constraint(cnst);
  free(jlmat);
  jlmat = NULL;

  aggregate_sparsecnst(ost,pi,ist->cnst,ist->k);
  ret = reduce_amortize(ost,pi);
  if(ret) {
    ret = 8;
    goto err;
  }

  polx_refresh(ost->cnst->b);
  polxvec_refresh(ost->cnst->phi,ost->n);
  return 0;

err:
  free_statement(ost);
  free_constraint(cnst);
  free(jlmat);
  return ret;
}

int principle_verify(const prncplstmnt *st, const witness *wt) {
  size_t i;
  int ret;
  uint64_t normsq = 0;
  polx *sx[wt->r];

  if(wt->r != st->r)
    return 1;
  for(i=0;i<wt->r;i++)
    if(wt->n[i] != st->n[i])
      return 2;
  for(i=0;i<wt->r;i++)
    normsq += polyvec_sprodz(wt->s[i],wt->s[i],wt->n[i]);
  if(normsq > st->betasq)
    return 3;

  *sx = NULL;
  for(i=0;i<st->k;i++) {
    ret = !sparsecnst_check(&st->cnst[i],sx,wt);
    if(ret) {
      ret = 10+i;
      goto end;
    }
  }

end:
  free(*sx);
  return ret;
}
