#include <stdlib.h>
#include <stddef.h>
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
#include "greyhound.h"

//#define STREAM_WITNESS

static int init_polcomctx(polcomctx *ctx, size_t len) {
  double varz;
  comparams *cpp = ctx->cpp;

  ctx->len = len;

  for(cpp->f=2;cpp->f<=8;cpp->f++) {
    cpp->b = (LOGQ+cpp->f/2)/cpp->f;
    for(cpp->kappa=1;cpp->kappa<=32;cpp->kappa++) {
      /* m*n = len
       * minimize 2*m*f + (kappa+1)*fu*n, we use f as an approximation for fu */
      ctx->m = round(sqrt(len*(cpp->kappa+1)/2.0));
      ctx->n = ceil((double)len/ctx->m);

      varz = exp2(2*cpp->b)/12*ctx->n*(TAU1+4*TAU2);
      cpp->bu = round(0.25*log2(12*varz));  // z (decomposed)
      //cpp->bu = round(0.5*cpp->b + 0.25*log2(ctx->n) + 0.25*log2(TAU1+4*TAU2)));
      cpp->fu = round((double)LOGQ/cpp->bu);

      ctx->normsq  = (exp2(2*cpp->bu)/12 + varz/exp2(2*cpp->bu))*ctx->m*cpp->f;
      ctx->normsq += (exp2(2*cpp->bu)*(cpp->fu - 1) + exp2(2*(LOGQ-(cpp->fu - 1)*cpp->bu)))/12*(cpp->kappa+1)*ctx->n;
      ctx->normsq *= N;

      if(sis_secure(cpp->kappa,6*T*SLACK*exp2(cpp->bu)*sqrt(ctx->normsq)))
        break;
    }
    for(cpp->kappa1=1;cpp->kappa1<=32;cpp->kappa1++)
      if(sis_secure(cpp->kappa1,2*SLACK*sqrt(ctx->normsq)))
        break;

    if(cpp->kappa <= 32 && cpp->kappa1 <= 32)
      break;
  }
  if(cpp->kappa > 32) {
    fprintf(stderr,"ERROR in init_polcomctx(): Cannot make inner commitments secure!\n");
    return 1;
  }
  if(cpp->kappa1 > 32) {
    fprintf(stderr,"ERROR in init_polcomctx(): Cannot make outer commitments secure!\n");
    return 2;
  }

  cpp->u1len = cpp->kappa*cpp->fu*ctx->n;
  cpp->u2len = cpp->fu*ctx->n;

  ctx->s  = NULL;
#ifdef STREAM_WITNESS
  ctx->sx = NULL;
#else
  ctx->sx = _aligned_alloc(64,ctx->m*cpp->f*ctx->n*sizeof(polx));
#endif
  ctx->t  = _aligned_alloc(64,cpp->u1len*sizeof(poly));
  ctx->u1 = _aligned_alloc(64,cpp->kappa1*sizeof(polz));
  return 0;
}

void free_polcomctx(polcomctx *ctx) {
  free(ctx->sx);
  free(ctx->t);
  //free(ctx->u1); // u1 free'd when freeing proof
  ctx->s = NULL;
  ctx->t = NULL;
}

static void init_polcomprf(polcomprf *pi, const polcomctx *ctx, uint32_t x, uint32_t y) {
  pi->len = ctx->len;
  pi->m = ctx->m;
  pi->n = ctx->n;
  pi->x = x;
  pi->y = y;
  *pi->cpp = *ctx->cpp;
  pi->u1 = ctx->u1;
  pi->u2 = _aligned_alloc(64,ctx->cpp->kappa1*sizeof(polz));
}

void free_polcomprf(polcomprf *pi) {
  free(pi->u1);
  free(pi->u2);
  pi->u1 = pi->u2 = NULL;
}

void print_polcomctx_pp(const polcomctx *ctx) {
  const comparams *cpp = ctx->cpp;

  printf("Polynomial commitment context:\n");
  printf("  Polynomial lengths: %zu\n",ctx->len);
  printf("  Lengths decomposition: %zu x %zu\n",ctx->m,ctx->n);
  printf("  Streaming polynomial: %s\n",(ctx->sx) ? "NO" : "YES");
  printf("  Commitment ranks: kappa = %zu; kappa1 = %zu\n",cpp->kappa,cpp->kappa1);
  printf("  Decomposition bases: b = %d; bu = %d\n",1 << cpp->b,1 << cpp->bu);
  printf("  Expansion factors: f = %zu; fu = %zu\n",cpp->f,cpp->fu);
  printf("  Predicted witness norm: %.2f\n",sqrt(ctx->normsq));
  printf("\n");
}

double print_polcomprf_pp(const polcomprf *pi) {
  double s;
  const comparams *cpp = pi->cpp;

  printf("Polynomial commitment evaluation proof:\n");
  printf("  Polynomial lengths: %zu\n",pi->len);
  printf("  Lengths decomposition: %zu x %zu\n",pi->m,pi->n);
  printf("  Evaluation point: %ld\n",pi->x);
  printf("  Commitment ranks: kappa = %zu; kappa1 = %zu\n",cpp->kappa,cpp->kappa1);
  printf("  Decomposition bases: b = %d; bu = %d\n",1 << cpp->b,1 << cpp->bu);
  printf("  Expansion factors: f = %zu; fu = %zu\n",cpp->f,cpp->fu);
  printf("  Witness norm: %.2f\n",sqrt(pi->normsq));
  s  = 2*cpp->kappa1*N*LOGQ;
  s /= 8192;
  printf("  Proof size: %.2f KB\n",s);
  printf("\n");
  return s;
}

int polcom_commit(polcomctx *ctx, const polz *s, size_t len) {
  int ret;

  ret = init_polcomctx(ctx,len);
  if(ret) return ret;  // commitments not secure (1/2)

  const size_t m = ctx->m;
  const size_t n = ctx->n;
  const comparams *cpp = ctx->cpp;

  size_t i,j;
  __attribute__((aligned(16)))
  uint8_t hashbuf[cpp->kappa1*N*QBYTES];
  polx t[cpp->kappa*cpp->fu];
  polx u[cpp->kappa1];

  init_comkey(MAX(m*cpp->f,n*extlen(cpp->kappa*cpp->fu,cpp->kappa1)));

  ctx->s = s;
#ifdef STREAM_WITNESS
  polx (*sx)[0] = _aligned_alloc(64,m*cpp->f*sizeof(polx));
#else
  polx (*sx)[m*cpp->f] = (polx (*)[m*cpp->f])ctx->sx;
#endif

  j = 0;
  for(i=0;i<n;i++) {
    /* inner commitments */
    polzvec_decompose_topolxvec(sx[i],&s[i*m],MIN(m,len-i*m),m,cpp->f,cpp->b);
    polxvec_mul_extension(t,comkey,sx[i],m*cpp->f,cpp->kappa,1);
    polxvec_decompose(&ctx->t[i*cpp->kappa*cpp->fu],t,cpp->kappa,cpp->fu,cpp->bu);

    /* outer commitment */
    polxvec_frompolyvec(t,&ctx->t[i*cpp->kappa*cpp->fu],cpp->kappa*cpp->fu);
    if(i == 0)
      j += polxvec_mul_extension(u,&comkey[j],t,cpp->kappa*cpp->fu,cpp->kappa1,1);
    else {
      j += polxvec_mul_extension(t,&comkey[j],t,cpp->kappa*cpp->fu,cpp->kappa1,1);
      polxvec_add(u,u,t,cpp->kappa1);
    }
  }

  if(!ctx->sx) free(sx);
  polzvec_frompolxvec(ctx->u1,u,cpp->kappa1);
  polzvec_bitpack(hashbuf,ctx->u1,cpp->kappa1);
  shake128(ctx->h,16,hashbuf,sizeof(hashbuf));
  return 0;
}

#if LOGQ <= 30
static int64_t mulmodq(int64_t a, int64_t b) {
  int64_t t,u;
  const int64_t q = ((int64_t)1 << LOGQ) - QOFF;
  const int64_t v = ((__int128)1 << (62+LOGQ))/q;
  //const __int128 off = (__int128)1 << (61+LOGQ);

  t = a*b;
  u = (__int128)t*v >> (62+LOGQ);
  u = t - u*q;
  return u;
}
#else
static int64_t mulmodq(int64_t a, int64_t b) {
  __int128 t;
  int64_t u;
  const int64_t q = ((int64_t)1 << LOGQ) - QOFF;
  const int64_t v = ((__int128)1 << (62+LOGQ))/q;
  //const __int128 off = (__int128)1 << (122-LOGQ);

  t = (__int128)a*b;
  u = t >> (2*LOGQ-61);
  u = (__int128)u*v >> (123-LOGQ);
  u = (int64_t)t - u*q;
  return u;
}
#endif

static int64_t expmodq(int64_t a, size_t k) {
  int64_t t;

  t = (k&1) ? a : 1;
  k >>= 1;
  while(k) {
    a = mulmodq(a,a);
    if(k&1) t = mulmodq(t,a);
    k >>= 1;
  }

  return t;
}

static void vec_scalemodq(int64_t r[N], const int64_t a[N], int64_t s) {
  size_t i;

  for(i=0;i<N;i++)
    r[i] = mulmodq(a[i],s);
/*
  __m512i f,g,h,k;
  const __m512i aa = _mm512_set1_epi64(a);
  const __m512i qinv = _mm512_set1_epi64(3210379595);
  const __m512i q = _mm512_set1_epi64((1LL << LOGQ) - QOFF);
  const __mmask16 mask = _cvtu32_mask16(0xAAAA);

  for(i=0;i<N/16;i++) {
    f = _mm512_load_si512((__m512i*)&vec[16*i]);
    f = _mm512_mul_epu32(f,aa);
    g = _mm512_srli_epi64(f,32);
    g = _mm512_mul_epu32(g,aa);
    h = _mm512_mul_epu32(f,qinv);
    k = _mm512_mul_epu32(g,qinv);
    h = _mm512_mul_epu32(h,q);
    k = _mm512_mul_epu32(k,q);
    f = _mm512_add_epi64(f,h);
    g = _mm512_add_epi64(g,k);
    f = _mm512_srli_epi64(f,32);
    f = _mm512_mask_blend_epi32(mask,f,g);
    _mm512_store_si512((__m512i*)&vec[16*i],f);
  }
*/
}

int64_t polzvec_eval(const polz *a, size_t len, int64_t x) {
  size_t i,j,k;
  int64_t t,y;

  y = 0;
  for(i=len;i>0;i--) {
    for(j=N;j>0;j--) {
      t = 0;
      for(k=0;k<L;k++)
        t += (int64_t)a[i-1].limbs[k].c[j-1] << 14*k;
      y  = mulmodq(y,x);
      y += t;
    }
  }

  y = mulmodq(y,1);
  return y;
}

void polcom_eval(witness *wt, polcomprf *pi, const polcomctx *ctx, int64_t x, int64_t y) {
  size_t i,j,k;
  const comparams *cpp = ctx->cpp;
  const size_t m = ctx->m;
  const size_t n = ctx->n;
  int64_t xx[m];  // xx[i] = x^(N*i)
  __attribute__((aligned(16)))
  uint8_t hashbuf[16+(2+cpp->kappa1*N)*QBYTES];

#ifdef STREAM_WITNESS
  polx (*sx)[0] = _aligned_alloc(64,m*cpp->f*sizeof(polx));
#else
  polx (*sx)[m*cpp->f] = (polx (*)[m*cpp->f])ctx->sx;
#endif

  /* public data hash */
  memcpy(hashbuf,ctx->h,16);
  hashbuf[16] = x >>  0;
  hashbuf[17] = x >>  8;
  hashbuf[18] = x >> 16;
  hashbuf[19] = x >> 24;
  hashbuf[20] = y >>  0;
  hashbuf[21] = y >>  8;
  hashbuf[22] = y >> 16;
  hashbuf[23] = y >> 24;

  /* init proof */
  init_polcomprf(pi,ctx,x,y);

  /* init witness */
  size_t wtn[4];
  wtn[0] = m*cpp->f;
  wtn[1] = wtn[0];
  wtn[2] = cpp->u1len;
  wtn[3] = cpp->u2len;
  init_witness_raw(wt,4,wtn);
  polyvec_copy(wt->s[2],ctx->t,wt->n[2]);

  /* powers of x */
  for(i=1;i<N;i*=2)
    x = mulmodq(x,x);
  xx[0] = 1;
  for(i=1;i<m;i++)
    xx[i] = mulmodq(xx[i-1],x);

  /* compute w */
  polx t[1],u[1];
  polz w[1];
  for(i=0;i<n;i++) {
/*
    polz_copy(w,&ctx->s[i*m]);
    for(j=1;j<m;j++)
      polz_scale_add(w,&ctx->s[i*m+j],xx[j]);
*/
    if(!ctx->sx)
      polzvec_decompose_topolxvec(sx[i],&ctx->s[i*m],MIN(m,ctx->len-i*m),m,cpp->f,cpp->b);
    for(k=0;k<cpp->f;k++) {
      polxvec_copy(t,&sx[i][k*m],1);
      for(j=1;j<m;j++)
        polx_scale_add(t,&sx[i][k*m+j],xx[j]);
      if(k == 0)
        polxvec_copy(u,t,1);
      else {
        polx_refresh(t);
        polx_scale_add(u,t,(int64_t)1 << k*cpp->b);
      }
    }
    polz_frompolx(w,u);
    polz_decompose(&wt->s[3][i],w,n,cpp->fu,cpp->bu);
  }

  /* commit to w */
  polx *tmp = _aligned_alloc(64,MAX(MAX(wt->n[3],cpp->kappa1),n+m*cpp->f)*sizeof(polx));
  polxvec_frompolyvec(tmp,wt->s[3],wt->n[3]);
  polxvec_mul_extension(tmp,comkey,tmp,wt->n[3],cpp->kappa1,1);
  polzvec_frompolxvec(pi->u2,tmp,cpp->kappa1);
  polzvec_bitpack(&hashbuf[24],pi->u2,cpp->kappa1);

  /* ammortize s */
  polx *c = &tmp[0];
  polx *z = &tmp[n];
  shake128(hashbuf,32,hashbuf,sizeof(hashbuf));
  polxvec_challenge(c,n,&hashbuf[16],0);
  if(!ctx->sx)
    polzvec_decompose_topolxvec(sx[0],&ctx->s[0],m,m,cpp->f,cpp->b);
  polxvec_polx_mul(z,&c[0],sx[0],m*cpp->f);
  for(i=1;i<n;i++) {
    if(!ctx->sx)
      polzvec_decompose_topolxvec(sx[i],&ctx->s[i*m],MIN(m,ctx->len-i*m),m,cpp->f,cpp->b);
    polxvec_polx_mul_add(z,&c[i],sx[i],m*cpp->f);
  }

  polxvec_decompose(wt->s[0],z,wt->n[0],2,cpp->bu);
  pi->normsq = 0;
  for(i=0;i<wt->r;i++) {
    wt->normsq[i] = polyvec_sprodz(wt->s[i],wt->s[i],wt->n[i]);
    pi->normsq += wt->normsq[i];
  }

  if(!ctx->sx) free(sx);
  free(tmp);
}

int polcom_reduce(prncplstmnt *st, const polcomprf *pi) {
  int ret;
  size_t i,j;
  const comparams *cpp = pi->cpp;
  const size_t m = pi->m;
  const size_t n = pi->n;
  int64_t x = pi->x;
  int64_t y = pi->y;
  int64_t s;
  int64_t xvec[N], xvec2[N];
  size_t stn[4], len;
  __attribute__((aligned(16)))
  uint8_t hashbuf[24+cpp->kappa1*N*QBYTES];
  polx *buf;

  if(!sis_secure(cpp->kappa,6*T*SLACK*exp2(cpp->bu)*sqrt(pi->normsq))) {
    fprintf(stderr,"ERROR in polcom_reduce(): Inner commitments not secure\n");
    return 1;
  }
  if(!sis_secure(cpp->kappa1,2*SLACK*sqrt(pi->normsq))) {
    fprintf(stderr,"ERROR in polcom_reduce(): Outer commitments not secure\n");
    return 2;
  }

  /* init principal statement */
  stn[0] = m*cpp->f;
  stn[1] = stn[0];
  stn[2] = cpp->u1len;
  stn[3] = cpp->u2len;
  ret = init_prncplstmnt_raw(st,4,stn,pi->normsq,5,0);
  if(ret)  // total witness norm too big
    return 3;

  init_comkey(MAX(m*cpp->f,n*extlen(cpp->kappa*cpp->fu,cpp->kappa1)));

  /* challenges */
  polx *c = _aligned_alloc(64,n*sizeof(polx));
  polzvec_bitpack(hashbuf,pi->u1,cpp->kappa1);
  shake128(hashbuf,16,hashbuf,cpp->kappa1*N*QBYTES);
  hashbuf[16] = x >>  0;
  hashbuf[17] = x >>  8;
  hashbuf[18] = x >> 16;
  hashbuf[19] = x >> 24;
  hashbuf[20] = y >>  0;
  hashbuf[21] = y >>  8;
  hashbuf[22] = y >> 16;
  hashbuf[23] = y >> 24;
  polzvec_bitpack(&hashbuf[24],pi->u2,cpp->kappa1);
  shake128(hashbuf,32,hashbuf,sizeof(hashbuf));
  memcpy(st->h,hashbuf,16);
  polxvec_challenge(c,n,&hashbuf[16],0);

  /* powers of x (sigmam1) */
  xvec[0] = 1;
  xvec[N-1] = -x;
  for(i=N-2;i>0;i--)
    xvec[i] = mulmodq(xvec[i+1],x);
  x = -mulmodq(xvec[1],x);  // -x^N

  /* outer commitments */
  sparsecnst *cnst = &st->cnst[0];
  init_sparsecnst_half(cnst,4,n,0,cpp->kappa1,0,0);
  polzvec_topolxvec(cnst->b,pi->u1,cpp->kappa1);
  len = cpp->fu*cpp->kappa;
  for(i=0;i<n;i++) {
    cnst->idx[i] = 2;
    cnst->off[i] = i*len;
    cnst->len[i] = len;
    cnst->mult[i] = 1;
    cnst->phi[i] = &comkey[i*extlen(len,cpp->kappa1)];
  }

  cnst = &st->cnst[1];
  init_sparsecnst_half(cnst,4,1,0,cpp->kappa1,0,0);
  polzvec_topolxvec(cnst->b,pi->u2,cpp->kappa1);
  cnst->idx[0] = 3;
  cnst->off[0] = 0;
  cnst->len[0] = st->n[3];
  cnst->mult[0] = 1;
  cnst->phi[0] = comkey;

  /* <al,z> = \sum_i=0^n-1 c_iw_i */
  cnst = &st->cnst[2];
  len = st->n[0]+st->n[1]+st->n[3];
  buf = init_sparsecnst_half(cnst,4,3,len,1,0,1);
  cnst->idx[0] = 0;
  cnst->off[0] = 0;
  cnst->len[0] = st->n[0];
  cnst->mult[0] = 1;
  cnst->phi[0] = buf;
  buf += st->n[0];
  for(i=0;i<cpp->f;i++) {
    s = (int64_t)1 << i*cpp->b;
    for(j=0;j<m;j++) {
      polx_monomial(&cnst->phi[0][i*m+j],s,0);
      s = mulmodq(s,x);
    }
  }

  cnst->idx[1] = 1;
  cnst->off[1] = 0;
  cnst->len[1] = st->n[1];
  cnst->mult[1] = 1;
  cnst->phi[1] = buf;
  buf += st->n[1];
  polxvec_scale(cnst->phi[1],cnst->phi[0],st->n[1],(int64_t)1 << cpp->bu);

  cnst->idx[2] = 3;
  cnst->off[2] = 0;
  cnst->len[2] = st->n[3];
  cnst->mult[2] = 1;
  cnst->phi[2] = buf;
  polxvec_neg(cnst->phi[2],c,n);
  for(i=1;i<cpp->fu;i++)
    polxvec_scale(&cnst->phi[2][i*n],&cnst->phi[2][0],n,(int64_t)1 << i*cpp->bu);

  /* inner commitments */
  /* Az = \sum_i=0^n-1 c_it_i */
  free(c);
  c = cnst->phi[2];
  cnst = &st->cnst[3];
  len = extlen(st->n[1],cpp->kappa) + st->n[2];
  buf = init_sparsecnst_half(cnst,4,3,len,cpp->kappa,0,1);
  cnst->idx[0] = 0;
  cnst->off[0] = 0;
  cnst->len[0] = st->n[0];
  cnst->mult[0] = 1;
  cnst->phi[0] = comkey;

  cnst->idx[1] = 1;
  cnst->off[1] = 0;
  cnst->len[1] = st->n[1];
  cnst->mult[1] = 1;
  cnst->phi[1] = buf;
  buf += extlen(st->n[1],cpp->kappa);
  polxvec_scale(cnst->phi[1],comkey,extlen(st->n[1],cpp->kappa),(int64_t)1 << cpp->bu);

  cnst->idx[2] = 2;
  cnst->off[2] = 0;
  cnst->len[2] = n*cpp->fu;
  cnst->mult[2] = cpp->kappa;
  cnst->phi[2] = buf;
  for(i=0;i<n;i++)
    for(j=0;j<cpp->fu;j++)
      polxvec_copy(&cnst->phi[2][(i*cpp->fu+j)*cpp->kappa],&c[n*j+i],1);

  /* <w,ar> = b */
  cnst = &st->cnst[4];
  buf = init_sparsecnst_half(cnst,4,1,st->n[3],0,0,0);
  cnst->idx[0] = 3;
  cnst->off[0] = 0;
  cnst->len[0] = st->n[3];
  cnst->mult[0] = 1;
  cnst->phi[0] = buf;
  polx_monomial(cnst->b,y,0);
  s = expmodq(x,m);
  for(i=0;i<cpp->fu;i++) {
    memcpy(xvec2,xvec,sizeof(xvec));
    for(j=0;j<n;j++) {
      polxvec_fromint64vec(&cnst->phi[0][i*n+j],1,1,xvec2);
      vec_scalemodq(xvec2,xvec2,s);
    }
    vec_scalemodq(xvec,xvec,(int64_t)1 << cpp->bu);
  }

  return 0;
}
