#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <stdio.h>
#include "malloc.h"
#include "fips202.h"
#include "data.h"
#include "poly.h"
#include "polx.h"
#include "polz.h"
#include "labrador.h"
#include "chihuahua.h"
#include "dachshund.h"

void init_smplstmnt_raw(smplstmnt *st, size_t r, const size_t n[r], const uint64_t betasq[r], size_t k) {
  size_t i;
  void *buf;

  buf = _malloc(r*sizeof(size_t) + k*sizeof(sparsecnst) + r*sizeof(uint64_t));
  st->n = (size_t*)buf;
  st->cnst = (sparsecnst*)&st->n[r];
  st->betasq = (uint64_t*)&st->cnst[k];

  st->r = r;
  st->k = k;
  for(i=0;i<r;i++) {
    st->n[i] = n[i];
    st->betasq[i] = betasq[i];
  }
  for(i=0;i<k;i++)
    st->cnst[i].idx = NULL;

  st->u0 = st->alpha = NULL;
  memset(st->h,0,16);
}

int set_smplstmnt_lincnst_raw(smplstmnt *st, size_t i, size_t nz, const size_t idx[nz],
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

void free_smplstmnt(smplstmnt *st) {
  size_t i;

  if(!st->n) return;
  for(i=0;i<st->k;i++)
    free_sparsecnst(&st->cnst[i]);
  free(st->n);  // one buffer for n,cnst,betasq
  free(st->u0);  // one buffer for u0,alpha
  st->n = NULL;
  st->u0 = NULL;
}

static void init_commitment(commitment *com, const comparams *cpp, const witness *wt) {
  size_t i;
  void *buf;

  com->r = wt->r;
  com->fu = cpp->fu;
  com->bu = cpp->bu;
  com->kappa = cpp->kappa;
  com->kappa1 = cpp->kappa1;
  buf = _malloc(com->r*sizeof(size_t));
  com->n = (size_t*)buf;
  for(i=0;i<com->r;i++)
    com->n[i] = wt->n[i];

  buf = _aligned_alloc(64,com->kappa1*sizeof(polz)+com->r*sizeof(polx));
  com->u = (polz*)buf;
  com->alpha = (polx*)&com->u[com->kappa1];
}

void free_commitment(commitment *com) {
  free(com->n);
  free(com->u);
  com->n = NULL;
  com->u = NULL;
}

void print_smplstmnt_pp(const smplstmnt *st) {
  size_t i;

  printf("Dachshund statement ");
  for(i=0;i<5;i++)
    printf("%02hhX",st->h[i]);
  printf(":\n");

  printf("  Witness multiplicity: %zu\n",st->r);
  printf("  Witness ranks: ");
  for(i=0;i<MIN(st->r,10);i++) {
    printf("%zu",st->n[i]);
    if(i<st->r-1) printf(", ");
  }
  printf("\n");

  printf("  Number of dot-product constraints: %zu\n",st->k);
  printf("  Norm constraints: ");
  for(i=0;i<MIN(st->r,10);i++) {
    printf("%.2f",sqrt(st->betasq[i]));
    if(i<st->r-1) printf(", ");
  }
  printf("\n\n");
}

static void expand_witness(witness *ewt, const smplstmnt *ist, const witness *iwt) {
  size_t i,k;
  const size_t r = iwt->r;
  void *buf;

  ewt->r = r+1;
  buf = _malloc(ewt->r*(sizeof(size_t) + sizeof(uint64_t) + sizeof(poly*)));
  ewt->n = (size_t*)buf;
  ewt->normsq = (uint64_t*)&ewt->n[ewt->r];
  ewt->s = (poly**)&ewt->normsq[ewt->r];
  k = 0;
  for(i=0;i<r;i++) {
    ewt->n[i] = iwt->n[i];
    ewt->normsq[i] = iwt->normsq[i];
    ewt->s[i] = iwt->s[i];
    if(ist->betasq[i]) k += 1;
  }

  buf = _aligned_alloc(64,k*sizeof(poly));
  ewt->s[r] = (poly*)buf;
  k = 0;
  for(i=0;i<r;i++) {
    if(!ist->betasq[i])
      continue;

    poly_binary_fromuint64(&ewt->s[r][k++],ist->betasq[i] - ewt->normsq[i]);
  }

  ewt->n[r] = k;
  ewt->normsq[r] = polyvec_sprodz(ewt->s[r],ewt->s[r],k);
/*
  uint64_t normsq = 0;
  for(i=0;i<r+1;i++)
    normsq += ewt->normsq[i];
  if(normsq >= (uint64_t)1 << LOGQ)
    printf("ERROR: Total witness norm too big to support exact l2-norm proofs\n");
*/
}

static void init_gadget(polx *gdgt) {
  size_t i;
  int64_t coeffs[N];

  for(i=0;i<30;i++)  // Maximum betasq = 2^30
    coeffs[i] = (int64_t)1 << i;
  while(i<N)
    coeffs[i++] = 0;
  polxvec_fromint64vec(gdgt,1,1,coeffs);
  polx_sigmam1(gdgt,gdgt);
}

static void simple_commit(statement *ost, witness *owt, proof *pi, commitment *com, polx sx[ost->r][ost->n],
                          const smplstmnt *ist, const witness *ewt)
{
  const size_t r = ewt->r;
  const size_t *n = ewt->n;
  const comparams *cpp = pi->cpp;

  const size_t s0 = 0;
  const size_t s1 = s0 + pi->nu[r-1];
  const size_t s2 = s1 + pi->nu[r];
  const size_t s3 = s2 + pi->nu[2*r];  // ost->r

  /* offsets in aux witness vector v */
  const size_t t = 0;
  const size_t t2 = t + s2*cpp->fu*cpp->kappa;
  const size_t g = t + s3*cpp->fu*cpp->kappa;

  size_t i,k,l;
  __attribute__((aligned(16)))
  uint8_t hashbuf[16+MAX(cpp->u1len,com->kappa1)*N*QBYTES];
  poly *v = owt->s[cpp->f];
  polx gdgt[1];
  polz tmp[1];

  init_gadget(gdgt);

  k = 0;
  for(i=0;i<r-1;i++) {
    polxvec_frompolyvec(&sx[s0][k],ewt->s[i],n[i]);
    if(ist->betasq[i])
      polxvec_sigmam1(&sx[s2][k],&sx[s0][k],n[i]);
    else
      polxvec_flip(&sx[s2][k],&sx[s0][k],n[i]);

    polxvec_sprod(&sx[s1][i],&sx[s0][k],&sx[s2][k],n[i]);
    k += n[i];
  }

  polxvec_frompolyvec(&sx[s0][k],ewt->s[r-1],n[r-1]);
  polxvec_flip(&sx[s2][k],&sx[s0][k],n[r-1]);
  polxvec_sprod(&sx[s1][r-1],&sx[s0][k],&sx[s2][k],n[r-1]);

  for(i=0;i<r;i++) {
    if(i < r-1 && ist->betasq[i])
      polx_mul_add(&sx[s1][i],gdgt,&sx[s0][k++]);
    polz_frompolx(tmp,&sx[s1][i]);
    polz_setcoeff_fromint64(tmp,0,0);
    polz_center(tmp);
    polz_decompose_topolx(&sx[s1][i],tmp,r,3,10);
  }

  polxvec_setzero(&sx[s0][k],(s1-s0)*ost->n - k);
  polxvec_setzero(&sx[s1][3*r],(s2-s1)*ost->n - 3*r);
  polxvec_setzero(&sx[s2][k],(s3-s2)*ost->n - k);

  polxvec_setzero(ost->u1,com->kappa1);
  l = commit_raw(ost->u1,&v[t],s2-s0,ost->n,&sx[s0],0,cpp,pi->tail);
  polzvec_frompolxvec(com->u,ost->u1,com->kappa1);
  memcpy(hashbuf,ost->h,16);
  polzvec_bitpack(&hashbuf[16],com->u,com->kappa1);
  shake128(hashbuf,32,hashbuf,16+com->kappa1*N*QBYTES);

  polxvec_challenge(com->alpha,r,&hashbuf[16],0);
  k = 0;
  for(i=0;i<r;i++) {
    polxvec_polx_mul(&sx[s2][k],&com->alpha[i],&sx[s2][k],n[i]);
    k += n[i];
  }

  l = commit_raw(ost->u1,&v[t2],s3-s2,ost->n,&sx[s2],l,cpp,pi->tail);
  qugarbage_raw(ost->u1,&v[g],ost->r,ost->n,sx,l,cpp,pi->tail);
  polzvec_frompolxvec(pi->u1,ost->u1,cpp->u1len);
  polzvec_bitpack(&hashbuf[16],pi->u1,cpp->u1len);
  shake128(ost->h,16,hashbuf,16+cpp->u1len*N*QBYTES);
}

/*
  for(i=0;i<iwt->r;i++) {
    b = round(log2(12.0*iwt->normsq[i]/(iwt->n[i]*N))/2);
    b = MAX(1,b);
  }
  t = ceil(log2(iwt->normsq[i])/b);
*/

static void reduce_simple_commit(statement *ost, const proof *pi, const commitment *com) {
  __attribute__((aligned(16)))
  uint8_t hashbuf[16+MAX(pi->cpp->u1len,com->kappa1)*N*QBYTES];

  polzvec_topolxvec(ost->u1,pi->u1,pi->cpp->u1len);

  memcpy(hashbuf,ost->h,16);
  polzvec_bitpack(&hashbuf[16],com->u,com->kappa1);
  shake128(hashbuf,32,hashbuf,16+com->kappa1*N*QBYTES);
  //FIXME: check alpha
  polzvec_bitpack(&hashbuf[16],pi->u1,pi->cpp->u1len);
  shake128(ost->h,16,hashbuf,16+pi->cpp->u1len*N*QBYTES);
}

static void simple_collaps(constraint *cnst, statement *ost,
                           const proof *pi, const commitment *com, const smplstmnt *ist)
{
  const size_t r = ist->r+1;
  const size_t s0 = 0;
  const size_t s1 = s0 + pi->nu[r-1];
  const size_t s2 = s1 + pi->nu[r];

  size_t i,j,k;
  int64_t s;
  uint8_t hashbuf[32+r*QBYTES];
  uint64_t nonce = 0;
  polx (*phi)[ost->n] = (polx(*)[ost->n])cnst->phi;
  polx ones[1];
  polx tmp[1];

  for(i=0;i<N;i++)
    ones->vec[K-1].vec->c[i] = 1;
  polx_frompoly(ones,&ones->vec[K-1]);

  shake128(hashbuf,sizeof(hashbuf),ost->h,16);
  memcpy(ost->h,hashbuf,16);

  /* conjuagates */
  k = 0;
  for(i=0;i<r;i++) {
    for(j=0;j<pi->n[i];j++) {
      polxvec_almostuniform(tmp,1,&hashbuf[16],nonce++);
      polx_add(&phi[s2][k],&phi[s2][k],tmp);
      polx_mul(tmp,tmp,&com->alpha[i]);
      polx_sigmam1(tmp,tmp);
      if(i < r-1 && ist->betasq[i])
        polx_sub(&phi[s0][k],&phi[s0][k],tmp);
      else {
        polx_add(&phi[s0][k],&phi[s0][k],tmp);
        polx_mul_add(cnst->b,ones,tmp);
      }
      k += 1;
    }
  }

  /* lifting polynomials */
  for(i=0;i<r;i++) {
    s = 0;
    for(j=0;j<QBYTES;j++)
      s |= (int64_t)hashbuf[32+i*QBYTES+j] << 8*j;
    s &= (1LL << LOGQ) - 1;

    for(j=0;j<3;j++) {
      polx_monomial(tmp,s << 10*j,0);
      polx_add(&phi[s1][j*r+i],&phi[s1][j*r+i],tmp);
    }
  }
}

static void simple_aggregate(statement *ost, const proof *pi, const commitment *com, const smplstmnt *ist) {
  const size_t r = ist->r+1;
  const size_t s0 = 0;
  const size_t s1 = s0 + pi->nu[r-1];
  const size_t s2 = s1 + pi->nu[r];

  size_t i,j,k;
  sparsemat *a = ost->cnst->a;
  polx (*phi)[ost->n] = (polx(*)[ost->n])ost->cnst->phi;
  polx *b = ost->cnst->b;
  polx one[1];
  polx gdgt[1];

  polx_monomial(one,1,0);
  init_gadget(gdgt);

  /* assumes that non-zero coeffs in a are in row-major order */
  j = 0;
  for(i=0;i<s1;i++) {
    while(j < a->len && (a->rows[j] != i || a->cols[j] != (s2+i)))
      j++;
    if(j < a->len)
      polx_add(&a->coeffs[j],&a->coeffs[j],one);
    else {
      polxvec_copy(&a->coeffs[j],one,1);
      a->rows[j] = i;
      a->cols[j] = s2+i;
      a->len += 1;
    }
  }

  k = 0;
  for(i=0;i<r-1;i++)
    k += pi->n[i];

  for(i=0;i<r;i++) {
    if(i < r-1 && ist->betasq[i]) {
      polx_scale_add(b,&com->alpha[i],ist->betasq[i]);
      polx_mul_add(&phi[s0][k++],gdgt,&com->alpha[i]);
    }

    polx_sub(&phi[s1][i],&phi[s1][i],&com->alpha[i]);
    for(j=1;j<3;j++)
      polx_scale_add(&phi[s1][j*r+i],&com->alpha[i],-(int64_t)(1 << 10*j));
  }
}

void simple_prove(statement *ost, witness *owt, proof *pi, commitment *com,
                  const smplstmnt *ist, const witness *iwt, int tail)
{
  size_t i;
  witness ewt[1];
  constraint cnst[1];

  expand_witness(ewt,ist,iwt);
  init_proof(pi,ewt,2,tail);
  init_commitment(com,pi->cpp,ewt);
  init_statement(ost,pi,ist->h);
  init_witness(owt,ost);
  printf("Predicted witness norm: %.2f\n\n",sqrt(pi->normsq));

  polx (*sx)[ost->n] = _aligned_alloc(64,ost->r*ost->n*sizeof(polx));
  simple_commit(ost,owt,pi,com,sx,ist,ewt);

  const size_t s1 = pi->nu[ist->r];
  uint8_t (*jlmat)[ost->n][256*N/8] = _aligned_alloc(64,s1*ost->n*256*N/8);
  project(ost,pi,jlmat,ewt);

  free(ewt->s[ewt->r-1]);
  free(ewt->n);

  init_constraint_raw(cnst,ost->r,ost->n,1,0);
  for(i=0;i<LIFTS;i++) {
    collaps_jlproj_raw(cnst,s1,ost->n,ost->h,pi->p,jlmat);
    polxvec_setzero(&cnst->phi[s1*ost->n],(ost->r-s1)*ost->n);
    collaps_sparsecnst(cnst,ost,pi,ist->cnst,ist->k);
    simple_collaps(cnst,ost,pi,com,ist);
    lift_aggregate_zqcnst(ost,pi,i,cnst,sx);
  }
  free(jlmat);
  free_constraint(cnst);

  simple_aggregate(ost,pi,com,ist);
  aggregate_sparsecnst(ost,pi,ist->cnst,ist->k);
  amortize(ost,owt,pi,sx);
  free(sx);

  polx_refresh(ost->cnst->b);
  polxvec_refresh(ost->cnst->phi,ost->n);
}

int simple_reduce(statement *ost, const proof *pi, const commitment *com, const smplstmnt *ist) {
  size_t i;
  int ret;
  uint64_t betasq = 0;
  uint8_t (*jlmat)[ost->n][256*N/8];
  constraint cnst[1];

  init_statement(ost,pi,ist->h);
  for(i=0;i<ist->r;i++) {
    betasq += ist->betasq[i];
    if(ist->betasq[i])
      betasq += 30;
  }

  const size_t s1 = pi->nu[ist->r];
  jlmat = _aligned_alloc(64,s1*ost->n*256*N/8);

  reduce_simple_commit(ost,pi,com);
  ret = reduce_project(ost,jlmat,pi,ist->r+1,betasq);
  if(ret) goto err;

  init_constraint(cnst,ost);
  for(i=0;i<LIFTS;i++) {
    collaps_jlproj_raw(cnst,s1,ost->n,ost->h,pi->p,jlmat);
    polxvec_setzero(&cnst->phi[s1*ost->n],(ost->r-s1)*ost->n);
    collaps_sparsecnst(cnst,ost,pi,ist->cnst,ist->k);
    simple_collaps(cnst,ost,pi,com,ist);
    ret = reduce_lift_aggregate_zqcnst(ost,pi,i,cnst);
    if(ret) {
      ret = 2+i;
      goto err;
    }
  }
  free_constraint(cnst);
  free(jlmat);
  jlmat = NULL;

  simple_aggregate(ost,pi,com,ist);
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

int simple_verify(const smplstmnt *st, const witness *wt) {
  int ret = 0;
  size_t i;
  polx *sx[wt->r];

  if(wt->r != st->r)
    return 1;
  for(i=0;i<wt->r;i++) {
    if(wt->n[i] != st->n[i])
      return 2;
    if(!st->betasq[i] && !polyvec_isbinary(wt->s[i],wt->n[i]))
      return 3;
    if(st->betasq[i] && wt->normsq[i] > st->betasq[i])
      return 4;
  }

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
