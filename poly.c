#include <stdint.h>
#include <string.h>
#include <x86intrin.h>
#include <math.h>
#include <complex.h>
#include "aesctr.h"
#include "fips202.h"
#include "data.h"
#include "poly.h"

static int16_t fpmul(int16_t a, int16_t b, const pdata *prime) {
  int32_t c;
  int16_t t;

  c = (int32_t)a*b;
  t = (int16_t)c*prime->pinv;
  t = (c - (int32_t)t*prime->p) >> 16;
  return t;
}

static int16_t fpred(int16_t a, const pdata *prime) {
  int16_t t;

  t  = ((int32_t)prime->v*a + (1<<26)) >> 27;
  t *= prime->p;
  return a - t;
}

void polyvec_setzero(poly *r, size_t len) {
  size_t i,j;
  const __m512i zero = _mm512_setzero_si512();

  for(i=0;i<len;i++)
    for(j=0;j<N/32;j++)
      _mm512_store_si512(&r[i].vec->v[j],zero);
}

int polyvec_isbinary(const poly *r, size_t len) {
  int64_t ret = 0;
  poly t[16];

  while(len >= 16) {
    polyvec_flip(t,r,16);
    ret |= polyvec_sprodz(t,r,16);

    r += 16;
    len -= 16;
  }
  polyvec_flip(t,r,len);
  ret |= polyvec_sprodz(t,r,len);

  ret = -ret >> 63;
  return ret;
}

void polyvec_fromint64vec(poly *r, size_t len, size_t deg, const int64_t a[len*deg*N]) {
  size_t i,j,k;

  for(i=0;i<len;i++)
    for(j=0;j<deg;j++)
      for(k=0;k<N;k++)
        r[i*deg+j].vec->c[k] = a[i*deg*N+k*deg+j];
}

void polyvec_copy(poly *r, const poly *a, size_t len) {
  size_t i,j;
  __m512i f;

  for(i=0;i<len;i++) {
    for(j=0;j<N/32;j++) {
      f = _mm512_load_si512(&a[i].vec->v[j]);
      _mm512_store_si512(&r[i].vec->v[j],f);
    }
  }
}

void poly_binary_fromuint64(poly *r, uint64_t a) {
  size_t i;

  for(i=0;i<N;i++) {
    r->vec->c[i] = a & 1;
    a >>= 1;
  }
}

static void poly_constant_ntt(poly *r, int16_t v, const pdata *prime) {
  size_t i;
  __m512i vv;

  v = fpmul(v,prime->montsq,prime);
  vv = _mm512_set1_epi16(v);
  for(i=0;i<N/32;i++)
    _mm512_store_si512(&r->vec->v[i],vv);
}

void poly_monomial_ntt(poly *r, int16_t v, int k, const pdata *prime) {
  if(k == 0) {
    poly_constant_ntt(r,v,prime);
    return;
  }

  polyvec_setzero(r,1);
  v = fpmul(v,prime->s,prime);
  r->vec->c[k] = v;
  poly_ntt(r,r,prime);
}

static inline size_t uniform_ref(int16_t *r, size_t len, const uint8_t *buf, size_t buflen, int16_t p) {
  size_t i,j;
  int s;
  int16_t t,u;
  const int bits = (p & 0x2000) ? 14 : 13;
  const int16_t mask = (1 << bits) - 1;

  i = j = 0;
  s = 0;
  while(i < len && j <= buflen-2) {
    t  = (s) ? (int16_t)buf[j-1] >> (8-s) : 0;
    t |= (int16_t)buf[j+0] << (s+0);
    t |= (int16_t)buf[j+1] << (s+8);
    t &= mask;

    s = s + 16 - bits;
    j += 2 - s/8;
    s %= 8;

    if(t < p) {
      u = ((p-1)/2 - t) >> 15;
      t -= u & p;
      r[i++] = t;
    }
  }

  return i;
}

static inline size_t uniform(int16_t *r, size_t len, const uint8_t *buf, size_t buflen, int16_t p) {
  size_t i,j;
  const int bits = (p & 0x2000) ? 14 : 13;
  __m512i f,g;
  __mmask32 store, center;
  const __m512i pp = _mm512_set1_epi16(p);
  const __m512i hp = _mm512_set1_epi16(p/2);
  const __m512i mask = _mm512_set1_epi16((1 << bits) - 1);
  __m512i permbidx,srlvwidx,sllvwidx;
  if(bits == 13) {
    permbidx = _mm512_set_epi8(51,50,49,48,48,47,46,45,44,43,43,42,41,40,40,39,
                               38,37,36,35,35,34,33,32,31,30,30,29,28,27,27,26,
                               25,24,23,22,22,21,20,19,18,17,17,16,15,14,14,13,
                               12,11,10, 9, 9, 8, 7, 6, 5, 4, 4, 3, 2, 1, 1, 0);
    srlvwidx = _mm512_broadcast_i64x2(_mm_set_epi16( 3, 6, 1, 4, 7, 2, 5, 0));
    sllvwidx = _mm512_broadcast_i64x2(_mm_set_epi16(13,10,13,12, 9,13,11,13));
  }
  else {
    permbidx = _mm512_set_epi8(55,54,53,52,51,50,50,49,48,47,46,45,44,43,43,42,
                               41,40,39,38,37,36,36,35,34,33,32,31,30,29,29,28,
                               27,26,25,24,23,22,22,21,20,19,18,17,16,15,15,14,
                               13,12,11,10, 9, 8, 8, 7, 6, 5, 4, 3, 2, 1, 1, 0);
    srlvwidx = _mm512_broadcast_i64x2(_mm_set_epi16( 2, 4, 6, 0, 2, 4, 6, 0));
    sllvwidx = _mm512_broadcast_i64x2(_mm_set_epi16(14,12,10,14,14,12,10,14));
  }

  i = j = 0;
  while(i+32 <= len && j+4*bits <= buflen) {
    f = _mm512_loadu_si512(&buf[j]);
    j += 4*bits;
    f = _mm512_permutexvar_epi8(permbidx,f);
    g = _mm512_alignr_epi8(f,f,2);
    f = _mm512_srlv_epi16(f,srlvwidx);
    g = _mm512_sllv_epi16(g,sllvwidx);
    f = _mm512_or_si512(f,g);
    f = _mm512_and_si512(f,mask);
    store = _mm512_cmp_epi16_mask(f,pp,1);
    center = _mm512_cmp_epi16_mask(hp,f,1);
    f = _mm512_mask_sub_epi16(f,center,f,pp);
    _mm512_mask_compressstoreu_epi16(&r[i],store,f);
    i += _popcnt32(_cvtmask32_u32(store));
  }

  return i+uniform_ref(r+i,len-i,buf+j,buflen-j,p);
}

void polyvec_uniform(poly *r, size_t len, const pdata *prime, const uint8_t seed[16], uint64_t nonce) {
  size_t k;
  int16_t *coeffs;
  const int16_t p = prime->p;
  const int bits = (p & 0x2000) ? 14 : 13;
  size_t nblocks = ((32*N*bits/8 << bits)/p + AES128CTR_BLOCKBYTES-1)/AES128CTR_BLOCKBYTES;
  __attribute__((aligned(64)))
  uint8_t buf[nblocks*AES128CTR_BLOCKBYTES];
  aes128ctr_ctx aesctx;

  aes128ctr_init(&aesctx,seed,nonce);

  coeffs = r->vec->c;
  len *= N;

  while(len >= 32*N) {
    aes128ctr_squeezeblocks(buf,nblocks,&aesctx);
    k = uniform(coeffs,32*N,buf,nblocks*AES128CTR_BLOCKBYTES,p);
    coeffs += k;
    len -= k;
  }

  while(len) {
    nblocks = ((len*bits/8 << bits)/p + AES128CTR_BLOCKBYTES-1)/AES128CTR_BLOCKBYTES;
    aes128ctr_squeezeblocks(buf,nblocks,&aesctx);
    k = uniform(coeffs,len,buf,nblocks*AES128CTR_BLOCKBYTES,p);
    coeffs += k;
    len -= k;
  }
}

/* len must be even */
static inline void ternary_ref(int16_t *r, size_t len, const uint8_t *buf) {
  size_t i;
  uint8_t t;
  const uint16_t lut = 0xA815;

  for(i=0;i<len/2;i++) {
    t = buf[i];
    r[2*i+0]  = (lut >> (t & 0xF)) & 0x3;
    r[2*i+0] -= 1;

    r[2*i+1]  = (lut >> (t >> 4)) & 0x3;
    r[2*i+1] -= 1;
  }
}

static inline void ternary(__m512i *r, size_t len, const uint8_t *buf) {
  size_t i;
  __m512i f,g;
  const __m512i one = _mm512_set1_epi16(1);
  const __m512i mask2 = _mm512_set1_epi16(0x3);
  const __m512i mask4 = _mm512_set1_epi16(0xF);
  const __m512i lut = _mm512_set1_epi16((int16_t)0xA815);

  for(i=0;i<len;i++) {
    f = _mm512_cvtepu8_epi32(_mm_load_si128((__m128i*)&buf[16*i]));
    g = _mm512_slli_epi32(f,12);
    f = _mm512_or_si512(f,g);
    f = _mm512_and_si512(f,mask4);
    f = _mm512_srlv_epi16(lut,f);
    f = _mm512_and_si512(f,mask2);
    f = _mm512_sub_epi16(f,one);
    _mm512_store_si512(&r[i],f);
  }
}

/* Samples ternary polynomial with probabilities 5/16,6/16,5/16 for -1,0,1, respectively (cbd(2) mod 3) */
void polyvec_ternary(poly *r, size_t len, const uint8_t seed[16], uint64_t nonce) {
  size_t nblocks;
  uint8_t *buf;
  aes128ctr_ctx aesctx;

  aes128ctr_init(&aesctx,seed,nonce);

  while(len >= 128) {
    nblocks = 64*N/AES128CTR_BLOCKBYTES;
    buf = (uint8_t*)r + 128*3*N/2;
    aes128ctr_squeezeblocks(buf,nblocks,&aesctx);
    ternary(r->vec->v,128*N/32,buf);
    len -= 128;
    r += 128;
  }

  if(len) {
    nblocks = (len*N/2 + AES128CTR_BLOCKBYTES-1)/AES128CTR_BLOCKBYTES;
    __attribute__((aligned(64)))
    uint8_t buf2[nblocks*AES128CTR_BLOCKBYTES];
    aes128ctr_squeezeblocks(buf2,nblocks,&aesctx);
    ternary(r->vec->v,len*N/32,buf2);
  }
}

static inline void quarternary(__m512i *r, size_t len, const uint8_t *buf) {
  size_t i;
  __m512i f,g;
  const __m512i two = _mm512_set1_epi16(2);
  const __m512i mask2 = _mm512_set1_epi16(0x3);

  for(i=0;i<len;i++) {
    f = _mm512_cvtepu8_epi64(_mm_loadl_epi64((__m128i*)&buf[8*i]));
    g = _mm512_slli_epi64(f,28);
    f = _mm512_or_si512(f,g);
    g = _mm512_slli_epi32(f,14);
    f = _mm512_or_si512(f,g);
    f = _mm512_and_si512(f,mask2);
    f = _mm512_sub_epi16(f,two);
    _mm512_store_si512(&r[i],f);
  }
}

void polyvec_quarternary(poly *r, size_t len, const uint8_t seed[16], uint64_t nonce) {
  size_t nblocks;
  uint8_t *buf;
  aes128ctr_ctx aesctx;

  aes128ctr_init(&aesctx,seed,nonce);

  while(len >= 256) {
    nblocks = 64*N/AES128CTR_BLOCKBYTES;
    buf = (uint8_t*)r + 256*7*N/4;
    aes128ctr_squeezeblocks(buf,nblocks,&aesctx);
    quarternary(r->vec->v,256*N/32,buf);
    len -= 256;
    r += 256;
  }

  if(len) {
    nblocks = (len*N/4 + AES128CTR_BLOCKBYTES-1)/AES128CTR_BLOCKBYTES;
    __attribute__((aligned(64)))
    uint8_t buf2[nblocks*AES128CTR_BLOCKBYTES];
    aes128ctr_squeezeblocks(buf2,nblocks,&aesctx);
    quarternary(r->vec->v,len*N/32,buf2);
  }
}

static inline size_t challenge(poly *c, size_t len, const uint8_t *buf, size_t buflen) {
  size_t i,j;
  int k,b;
  uint64_t signs;

  i = j = 0;
  while(i < len && j <= buflen-(TAU1+TAU2+(TAU1+TAU2+7)/8)) {
    signs = 0;
    for(k=0;k<(TAU1+TAU2+7)/8;k++)
      signs |= (uint64_t)buf[j++] << 8*k;

    polyvec_setzero(&c[i],1);

    k = N-TAU1-TAU2;
    while(k < N && j < buflen) {
      b = buf[j++] & (N-1);
      if(b <= k) {
        c[i].vec->c[k] = c[i].vec->c[b];
        c[i].vec->c[b] = (k < N-TAU2) ? 1 : 2;
        c[i].vec->c[b] -= (signs & 1) * (2*c[i].vec->c[b]);
        signs >>= 1;
        k += 1;
      }
    }

    if(k == N && poly_opnorm(&c[i]) <= T)
      i += 1;
  }

  return i;
}

void polyvec_challenge(poly *c, size_t len, const uint8_t seed[16], uint64_t nonce) {
  size_t k;
  __attribute__((aligned(64)))
  uint8_t buf[17*SHAKE128_RATE];
  shake128incctx shakectx;

  for(k=0;k<8;k++)
    buf[k] = nonce >> 8*k;

  shake128_inc_init(&shakectx);
  shake128_inc_absorb(&shakectx,seed,16);
  shake128_inc_absorb(&shakectx,buf,8);
  shake128_inc_finalize(&shakectx);

  while(len >= 10) {
    shake128_inc_squeezeblocks(buf,17,&shakectx);
    k = challenge(c,10,buf,17*SHAKE128_RATE);
    len -= k;
    c += k;
  }

  while(len) {
    k = (len*17+9)/10;
    shake128_inc_squeezeblocks(buf,k,&shakectx);
    k = challenge(c,len,buf,k*SHAKE128_RATE);
    len -= k;
    c += k;
  }
}

int64_t polyvec_sprodz_ref(const poly *a, const poly *b, size_t len) {
  size_t i,j;
  int64_t t=0;

  for(i=0;i<len;i++)
    for(j=0;j<N;j++)
      t += (int64_t)a[i].vec->c[j]*b[i].vec->c[j];

  return t;
}

int64_t polyvec_sprodz(const poly *a, const poly *b, size_t len) {
  size_t i,j;
  __m512i f,g,h,acc;
  __m256i t;
  __m128i u;

  acc = _mm512_setzero_si512();
  for(i=0;i<len;i++) {
    h = _mm512_setzero_si512();
    for(j=0;j<N/32;j++) {
      f = _mm512_load_si512(&a[i].vec->v[j]);
      g = _mm512_load_si512(&b[i].vec->v[j]);
      h = _mm512_dpwssd_epi32(h,f,g);
    }
    f = (__m512i)_mm512_moveldup_ps((__m512)h);
    g = _mm512_srai_epi64(h,32);
    f = _mm512_srai_epi64(f,32);
    acc = _mm512_add_epi64(acc,g);
    acc = _mm512_add_epi64(acc,f);
  }

  t = _mm256_add_epi64(_mm512_castsi512_si256(acc),_mm512_extracti64x4_epi64(acc,1));
  u = _mm_add_epi64(_mm256_castsi256_si128(t),_mm256_extracti64x2_epi64(t,1));
  u = _mm_add_epi64(u,_mm_unpackhi_epi64(u,u));

  return _mm_extract_epi64(u,0);
}

double polyvec_norm(const poly *a, size_t len) {
  return sqrt((double)polyvec_sprodz(a,a,len));
}

void poly_reduce(poly *r, const pdata *prime) {
  size_t i;
  __m512i f,g;
  const __m512i p = _mm512_set1_epi16(prime->p);
  const __m512i v = _mm512_set1_epi16(prime->v);
  const __m512i shift = _mm512_set1_epi16(1 << (16+15-27));

  for(i=0;i<N/32;i++) {
    f = _mm512_load_si512(&r->vec->v[i]);
    g = _mm512_mulhi_epi16(f,v);
    g = _mm512_mulhrs_epi16(g,shift);
    g = _mm512_mullo_epi16(g,p);
    f = _mm512_sub_epi16(f,g);
    _mm512_store_si512(&r->vec->v[i],f);
  }
}

void polyvec_reduce(poly *r, size_t len, size_t stride, const pdata *prime) {
  size_t i;

  for(i=0;i<len;i++)
    poly_reduce(&r[stride*i],prime);
}

void poly_center(poly *r, const pdata *prime) {
  size_t i;
  __m512i f;
  __mmask32 mask;
  const __m512i p = _mm512_set1_epi16(prime->p);
  const __m512i hp = _mm512_srli_epi16(p,1);

  for(i=0;i<N/32;i++) {
    f = _mm512_load_si512(&r->vec->v[i]);
    mask = _mm512_cmp_epi16_mask(hp,f,1);
    f = _mm512_mask_sub_epi16(f,mask,f,p);
    _mm512_store_si512(&r->vec->v[i],f);
  }
}

void polyvec_center(poly *r, size_t len, size_t stride, const pdata *prime) {
  size_t i;

  for(i=0;i<len;i++)
    poly_center(&r[stride*i],prime);
}

void poly_csubp(poly *r, const pdata *prime) {
  size_t i;
  __m512i f;
  __mmask32 mask;
  const __m512i p = _mm512_set1_epi16(prime->p);

  for(i=0;i<N/32;i++) {
    f = _mm512_load_si512(&r->vec->v[i]);
    mask = _mm512_cmp_epi16_mask(p,f,2);
    f = _mm512_mask_sub_epi16(f,mask,f,p);
    _mm512_store_si512(&r->vec->v[i],f);
  }
}

void polyvec_csubp(poly *r, size_t len, size_t stride, const pdata *prime) {
  size_t i;

  for(i=0;i<len;i++)
    poly_csubp(&r[stride*i],prime);
}

void poly_caddp(poly *r, const pdata *prime) {
  size_t i;
  __m512i f;
  __mmask32 mask;
  const __m512i p = _mm512_set1_epi16(prime->p);

  for(i=0;i<N/32;i++) {
    f = _mm512_load_si512(&r->vec->v[i]);
    mask = _mm512_movepi16_mask(f);
    f = _mm512_mask_add_epi16(f,mask,f,p);
    _mm512_store_si512(&r->vec->v[i],f);
  }
}

void polyvec_caddp(poly *r, size_t len, size_t stride, const pdata *prime) {
  size_t i;

  for(i=0;i<len;i++)
    poly_caddp(&r[stride*i],prime);
}

void poly_quot_add(poly *r, const poly *a, const pdata *prime) {
  size_t i;
  __m512i f,g;
  const __m512i v = _mm512_set1_epi16(prime->v);

  for(i=0;i<N/32;i++) {
    f = _mm512_load_si512(&a->vec->v[i]);
    g = _mm512_load_si512(&r->vec->v[i]);
    f = _mm512_mulhi_epi16(f,v);
    f = _mm512_add_epi16(f,g);
    _mm512_store_si512(&r->vec->v[i],f);
  }
}

void polyvec_quot_add(poly *r, const poly *a, size_t len, size_t stride, const pdata *prime) {
  size_t i;

  for(i=0;i<len;i++)
    poly_quot_add(&r[i],&a[stride*i],prime);
}

void poly_neg(poly *r, const poly *a) {
  size_t i;
  __m512i f;
  const __m512i zero = _mm512_setzero_si512();

  for(i=0;i<N/32;i++) {
    f = _mm512_load_si512(&a->vec->v[i]);
    f = _mm512_sub_epi16(zero,f);
    _mm512_store_si512(&r->vec->v[i],f);
  }
}

void polyvec_neg(poly *r, const poly *a, size_t len) {
  size_t i;

  for(i=0;i<len;i++)
    poly_neg(&r[i],&a[i]);
}

void poly_add(poly *r, const poly *a, const poly *b) {
  size_t i;
  __m512i f,g;

  for(i=0;i<N/32;i++) {
    f = _mm512_load_si512(&a->vec->v[i]);
    g = _mm512_load_si512(&b->vec->v[i]);
    f = _mm512_add_epi16(f,g);
    _mm512_store_si512(&r->vec->v[i],f);
  }
}

void polyvec_add(poly *r, const poly *a, const poly *b, size_t len) {
  size_t i;

  for(i=0;i<len;i++)
    poly_add(&r[i],&a[i],&b[i]);
}

void poly_sub(poly *r, const poly *a, const poly *b) {
  size_t i;
  __m512i f,g;

  for(i=0;i<N/32;i++) {
    f = _mm512_load_si512(&a->vec->v[i]);
    g = _mm512_load_si512(&b->vec->v[i]);
    f = _mm512_sub_epi16(f,g);
    _mm512_store_si512(&r->vec->v[i],f);
  }
}

void polyvec_sub(poly *r, const poly *a, const poly *b, size_t len) {
  size_t i;

  for(i=0;i<len;i++)
    poly_sub(&r[i],&a[i],&b[i]);
}

void poly_ntt_ref(poly *r, const pdata *prime) {
  int len, start, j, k;
  int16_t t;

  for(j=0;j<N;j++)
    r->vec->c[j] = fpmul(r->vec->c[j],prime->twist64[j],prime);

  for(len=N/2;len>=1;len>>=1) {
    for(start=0;start<N;start=j+len) {
      k = 0;
      for(j=start;j<start+len;j++) {
        t = r->vec->c[len+j];
        r->vec->c[len+j] = fpmul(r->vec->c[j] - t,prime->twist64[k],prime);
        r->vec->c[j] = fpred(r->vec->c[j] + t,prime);
        k += N/len;
      }
    }
  }
}

/*
void poly_ntt_old(poly * restrict r,const pdata *prime) {
  int len, start, j, k;
  int16_t t, zeta;

  k = 1;
  for(len=N/2;len>=1;len>>=1) {
    for(start=0;start<N;start=j+len) {
      zeta = prime->zetas[k++];
      for(j=start;j<start+len;j++) {
        t = fpmul(r->coeffs[j+len],zeta,prime);
        r->coeffs[j+len] = fpred(r->coeffs[j] - t,prime);
        r->coeffs[j] = fpred(r->coeffs[j] + t,prime);
      }
    }
  }
}

void poly_invntt_old(poly * restrict r,const pdata *prime) {
  int start, len, j, k;
  int16_t t, zeta;

  k = N-1;
  for(len=1;len<=N/2;len<<=1) {
    for(start=0;start<N;start=j+len) {
      zeta = prime->zetas[k--];
      for(j=start;j<start+len;j++) {
        t = r->coeffs[j];
        r->coeffs[j] = fpred(t + r->coeffs[j+len],prime);
        r->coeffs[j+len] = fpmul(r->coeffs[j+len] - t,zeta,prime);
      }
    }
  }

  poly_scale(r,r,prime,prime->f);
}
*/

void polyvec_ntt(poly *r, const poly *a, size_t len, size_t stride, const pdata *prime) {
  size_t i;

  //TODO: NTTx8

  for(i=0;i<len;i++)
    poly_ntt(&r[stride*i],&a[stride*i],prime);
}

void polyvec_invntt(poly *r, const poly *a, size_t len, size_t stride, const pdata *prime) {
  size_t i;

  for(i=0;i<len;i++)
    poly_invntt(&r[stride*i],&a[stride*i],prime);
}

void poly_pointwise(poly *r, const poly *a, const poly *b, const pdata *prime) {
  size_t i;
  __m512i f,g,h;
  const __m512i p = _mm512_set1_epi16(prime->p);
  const __m512i pinv = _mm512_set1_epi16(prime->pinv);

  for(i=0;i<N/32;i++) {
    f = _mm512_load_si512(&a->vec->v[i]);
    g = _mm512_load_si512(&b->vec->v[i]);
    h = _mm512_mullo_epi16(f,g);
    f = _mm512_mulhi_epi16(f,g);
    g = _mm512_mullo_epi16(h,pinv);
    g = _mm512_mulhi_epi16(g,p);
    f = _mm512_sub_epi16(f,g);
    _mm512_store_si512(&r->vec->v[i],f);
  }
}

void polyvec_pointwise(poly *r, const poly *a, const poly *b, size_t len, size_t stride, const pdata *prime) {
  size_t i;

  for(i=0;i<len;i++)
    poly_pointwise(&r[stride*i],&a[stride*i],&b[stride*i],prime);
}

void polyvec_poly_pointwise(poly *r, const poly *a, const poly *b, size_t len, size_t stride, const pdata *prime) {
  size_t i,j;
  __m512i f,g;
  __m512i al[N/32], ah[N/32];
  const __m512i p = _mm512_set1_epi16(prime->p);
  const __m512i pinv = _mm512_set1_epi16(prime->pinv);

  for(i=0;i<N/32;i++) {
    ah[i] = _mm512_load_si512(&a->vec->v[i]);
    al[i] = _mm512_mullo_epi16(ah[i],pinv);
  }

  for(i=0;i<len;i++) {
    for(j=0;j<N/32;j++) {
      f = _mm512_load_si512(&b[stride*i].vec->v[j]);
      g = _mm512_mullo_epi16(f,al[j]);
      f = _mm512_mulhi_epi16(f,ah[j]);
      g = _mm512_mulhi_epi16(g,p);
      f = _mm512_sub_epi16(f,g);
      _mm512_store_si512(&r[stride*i].vec->v[j],f);
    }
  }
}

void poly_pointwise_add(poly *r, const poly *a, const poly *b, const pdata *prime) {
  poly t;

  poly_pointwise(&t,a,b,prime);
  poly_add(r,r,&t);
}

void polyvec_pointwise_add(poly *r, const poly *a, const poly *b, size_t len, size_t stride, const pdata *prime) {
  size_t i;

  for(i=0;i<len;i++)
    poly_pointwise_add(&r[stride*i],&a[stride*i],&b[stride*i],prime);
}

// TODO: precompute a*qinv
void polyvec_poly_pointwise_add(poly *r, const poly *a, const poly *b, size_t len, size_t stride, const pdata *prime) {
  size_t i;

  for(i=0;i<len;i++)
    poly_pointwise_add(&r[stride*i],a,&b[stride*i],prime);
}

void polyvec_sprod_pointwise(poly *r, const poly *a, const poly *b, size_t len, size_t stride, const pdata *prime) {
  size_t i;
  const int extrared = (prime->p & 0x2000) ? 1 : 0;

  if(!len) {
    polyvec_setzero(r,1);
    return;
  }

  poly_pointwise(r,&a[0],&b[0],prime);
  for(i=1;i<len-1;i+=2) {
    poly_pointwise_add(r,&a[stride*(i+0)],&b[stride*(i+0)],prime);
    if(extrared) poly_reduce(r,prime);
    poly_pointwise_add(r,&a[stride*(i+1)],&b[stride*(i+1)],prime);
    poly_reduce(r,prime);
  }
  if(i<len) {
    poly_pointwise_add(r,&a[stride*i],&b[stride*i],prime);
    poly_reduce(r,prime);
  }
}

void polyvec_sprod_pointwise_add(poly *r, const poly *a, const poly *b, size_t len, size_t stride, const pdata *prime) {
  size_t i;
  const int extrared = (prime->p & 0x2000) ? 1 : 0;

  if(!len) return;

  for(i=0;i<len-1;i+=2) {
    poly_pointwise_add(r,&a[stride*(i+0)],&b[stride*(i+0)],prime);
    if(extrared) poly_reduce(r,prime);
    poly_pointwise_add(r,&a[stride*(i+1)],&b[stride*(i+1)],prime);
    poly_reduce(r,prime);
  }
  if(i<len) {
    poly_pointwise_add(r,&a[stride*i],&b[stride*i],prime);
    poly_reduce(r,prime);
  }
}

static size_t next2power(size_t a) {
  a -= 1;
  a |= a >>  1;
  a |= a >>  2;
  a |= a >>  4;
  a |= a >>  8;
  a |= a >> 16;
  a |= a >> 32;
  a += 1;
  return a;
}

size_t extlen(size_t len, size_t deg) {
  size_t mask;

  if(deg == 1)
    return len;

  mask = next2power(deg) - 1;
  return (len + mask) & ~mask;
}

static size_t bitrev(size_t k, size_t n) {
  const size_t t[32] = { 0, 16,  8, 24,  4, 20, 12, 28,
                         2, 18, 10, 26,  6, 22, 14, 30,
                         1, 17,  9, 25,  5, 21, 13, 29,
                         3, 19, 11, 27,  7, 23, 15, 31};

  return t[k] >> (5-n);
}

size_t polyvec_pointwise_extension(poly *c, const poly *a, const poly *b, size_t len, size_t stride, size_t deg,
                                   const pdata *prime)
{
  size_t i,j,k;
  size_t deg2;
  const int extrared = (prime->p & 0x2000) ? 1 : 0;
  poly tu[deg], tl[deg], tmp[stride*deg];

  /* shortcut */
  if(deg == 1) {
    polyvec_sprod_pointwise(c,a,b,len,stride,prime);
    return len;
  }

  deg2 = next2power(deg);
  polyvec_setzero(tu,deg);
  polyvec_setzero(tl,deg);

  k = 0;
  while(len) {
    for(i=0;i<MIN(deg,len);i++) {  // columns
      polyvec_poly_pointwise(&tmp[stride*0],&b[stride*i],&a[stride*(deg2-i)],i,stride,prime);
      polyvec_poly_pointwise(&tmp[stride*i],&b[stride*i],&a[stride*0],deg-i,stride,prime);
      for(j=0;j<i;j++)
        poly_add(&tu[j],&tu[j],&tmp[stride*j]);
      for(j=i;j<deg;j++)
        poly_add(&tl[j],&tl[j],&tmp[stride*j]);
      if(1 || extrared || i%2) { // FIXME
        polyvec_reduce(&tu[0],i,1,prime);
        polyvec_reduce(&tl[i],deg-i,1,prime);
      }
    }
    for(i=deg;i<MIN(deg2,len);i++) {
      polyvec_poly_pointwise(&tmp[stride*0],&b[stride*i],&a[stride*(deg2-i)],deg,stride,prime);
      for(j=0;j<deg;j++)
        poly_add(&tu[j],&tu[j],&tmp[stride*j]);
      if(1 || extrared || i%2) // FIXME
        polyvec_reduce(&tu[0],deg,1,prime);
    }

    a += stride*deg2;
    b += stride*deg2;
    k += deg2;
    len -= MIN(deg2,len);
  }

  poly x[1];
  poly_monomial_ntt(x,1,1,prime);
  polyvec_poly_pointwise(tu,x,tu,deg,1,prime);
  //polyvec_reduce(tl,deg,1,prime);
  for(i=0;i<deg;i++)
    poly_add(&c[stride*i],&tu[i],&tl[i]);
  polyvec_reduce(c,deg,stride,prime);

  return k;
}

size_t polyvec_collaps_add_extension(poly *c, const poly *a, const poly *b, size_t len, size_t stride, size_t deg,
                                     const pdata *prime)
{
  size_t i,k;
  size_t deg2;
  const poly *tl = a;
  poly tu[stride*deg];

  /* shortcut */
  if(deg == 1) {
    polyvec_poly_pointwise_add(c,a,b,len,stride,prime);
    return len;
  }

  deg2 = next2power(deg);
  poly x[1];
  poly_monomial_ntt(x,1,1,prime);
  polyvec_poly_pointwise(tu,x,tl,deg,stride,prime);

  k = 0;
  while(len) {
    for(i=0;i<MIN(deg,len);i++) {
      polyvec_sprod_pointwise_add(&c[stride*i],&tu[stride*0],&b[stride*(deg2-i)],i,stride,prime);
      polyvec_sprod_pointwise_add(&c[stride*i],&tl[stride*i],&b[stride*0],deg-i,stride,prime);
    }
    for(i=deg;i<MIN(deg2,len);i++)
      polyvec_sprod_pointwise_add(&c[stride*i],&tu[stride*0],&b[stride*(deg2-i)],deg,stride,prime);

    c += stride*deg2;
    b += stride*deg2;
    k += deg2;
    len -= MIN(deg2,len);
  }

  return k;
}

void poly_scale(poly *r, const poly *a, int16_t s, const pdata *prime) {
  size_t i;
  __m512i f,g;
  const __m512i l = _mm512_set1_epi16(s*prime->pinv);
  const __m512i h = _mm512_set1_epi16(s);
  const __m512i p = _mm512_set1_epi16(prime->p);

  for(i=0;i<N/32;i++) {
    f = _mm512_load_si512(&a->vec->v[i]);
    g = _mm512_mullo_epi16(f,l);
    f = _mm512_mulhi_epi16(f,h);
    g = _mm512_mulhi_epi16(g,p);
    f = _mm512_sub_epi16(f,g);
    _mm512_store_si512(&r->vec->v[i],f);
  }
}

void polyvec_scale(poly *r, const poly *a, size_t len, size_t stride, int16_t s, const pdata *prime) {
  size_t i;

  for(i=0;i<len;i++)
    poly_scale(&r[stride*i],&a[stride*i],s,prime);
}

void polyvec_scale_widening(poly *r, const poly *a, size_t len, size_t stride, int16_t s, const pdata *prime) {
  size_t i;

  for(i=0;i<len;i++)
    poly_scale(&r[stride*i],&a[i],s,prime);
}

void poly_scale_add(poly *r, const poly *a, int16_t s, const pdata *prime) {
  poly t;

  poly_scale(&t,a,s,prime);
  poly_add(r,r,&t);
}

void polyvec_scale_add(poly *r, const poly *a, size_t len, size_t stride, int16_t s, const pdata *prime) {
  size_t i;

  for(i=0;i<len;i++)
    poly_scale_add(&r[stride*i],&a[stride*i],s,prime);
}

static const double complex czetas[N/2] = {
   0.0                     + 1.0                    *I,
   0.70710678118654752440  + 0.70710678118654752440 *I,
   0.92387953251128675613  + 0.38268343236508977173 *I,
  -0.38268343236508977173  + 0.92387953251128675613 *I,
   0.98078528040323044913  + 0.19509032201612826785 *I,
  -0.19509032201612826785  + 0.98078528040323044913 *I,
   0.55557023301960222474  + 0.83146961230254523708 *I,
  -0.83146961230254523708  + 0.55557023301960222474 *I,
   0.99518472667219688624  + 0.098017140329560601994*I,
  -0.098017140329560601994 + 0.99518472667219688624 *I,
   0.63439328416364549822  + 0.77301045336273696081 *I,
  -0.77301045336273696081  + 0.63439328416364549822 *I,
   0.88192126434835502971  + 0.47139673682599764856 *I,
  -0.47139673682599764856  + 0.88192126434835502971 *I,
   0.29028467725446236764  + 0.95694033573220886494 *I,
  -0.95694033573220886494  + 0.29028467725446236764 *I,
   0.99879545620517239271  + 0.049067674327418014255*I,
  -0.049067674327418014255 + 0.99879545620517239271 *I,
   0.67155895484701840063  + 0.74095112535495909118 *I,
  -0.74095112535495909118  + 0.67155895484701840063 *I,
   0.90398929312344333159  + 0.42755509343028209432 *I,
  -0.42755509343028209432  + 0.90398929312344333159 *I,
   0.33688985339222005069  + 0.94154406518302077841 *I,
  -0.94154406518302077841  + 0.33688985339222005069 *I,
   0.97003125319454399260  + 0.24298017990326388995 *I,
  -0.24298017990326388995  + 0.97003125319454399260 *I,
   0.51410274419322172659  + 0.85772861000027206990 *I,
  -0.85772861000027206990  + 0.51410274419322172659 *I,
   0.80320753148064490981  + 0.59569930449243334347 *I,
  -0.59569930449243334347  + 0.80320753148064490981 *I,
   0.14673047445536175166  + 0.98917650996478097345 *I,
  -0.98917650996478097345  + 0.14673047445536175166 *I,
};

void poly_fft(double complex r[N/2], const poly *a) {
  size_t len, start, j, k;
  double complex t;

  for(j=0;j<N/2;j++)
    r[j] = a->vec->c[j] + a->vec->c[j]*I;

  k = 1;
  for(len=N/4;len>=1;len>>=1) {
    for(start=0;start<N/2;start=j+len) {
      for(j=start;j<start+len;j++) {
        t = r[j+len]*czetas[k];
        r[j+len] = r[j] - t;
        r[j] = r[j] + t;
      }
      k += 1;
    }
  }
}

void poly_invfft(poly *r, double complex a[N/2]) {
  size_t start, len, j, k;
  double complex u;

  k = N/4;
  for(len=1;len<=N/4;len<<=1) {
    for(start=0;start<N/2;start=j+len) {
      for(j=start;j<start+len;j++) {
        u = a[j] - a[j+len];
        a[j] = a[j] + a[j+len];
        a[j+len] = u*conj(czetas[k]);
      }
      k += 1;
    }
    k /= 4;
  }

  for(j=0;j<N/2;j++) {
    a[j] = ldexp(a[j],-5);
    r->vec->c[j] = round(creal(a[j]));
    r->vec->c[j+N/2] = round(cimag(a[j]));
  }
}

double poly_opnorm(const poly *a) {
  size_t i;
  double complex vec[N/2];
  double t,r = 0;

  poly_fft(vec,a);
  for(i=0;i<N/2;i++) {
    t = cabs(vec[i]);
    if(t > r) r = t;
  }

  return r;
}

void poly_sigmam1(poly *r, const poly *a) {
  size_t i;
  __m512i f,g;
  const __m512i permwidx = _mm512_set_epi16( 0, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12,13,14,15,
                                            16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31);
  const __m512i zero = _mm512_setzero_si512();

  r->vec->c[0] = a->vec->c[0];
  for(i=0;i<N/64;i++) {
    f = _mm512_loadu_si512(&a->vec->c[32*i+1]);
    g = _mm512_load_si512(&a->vec->v[N/32-1-i]);
    f = _mm512_permutexvar_epi16(permwidx,f);
    g = _mm512_permutexvar_epi16(permwidx,g);
    f = _mm512_sub_epi16(zero,f);
    g = _mm512_sub_epi16(zero,g);
    _mm512_storeu_si512(&r->vec->c[32*i+1],g);
    _mm512_store_si512(&r->vec->v[N/32-1-i],f);
  }
}

void polyvec_sigmam1(poly *r, const poly *a, size_t len) {
  size_t i;

  for(i=0;i<len;i++)
    poly_sigmam1(&r[i],&a[i]);
}

void poly_sigmam1_ntt(poly *r, const poly *a) {
  size_t i;
  __m512i f,g;
  const __m512i permwidx = _mm512_set_epi16( 0, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12,13,14,15,
                                            16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31);

  for(i=0;i<N/64;i++) {
    f = _mm512_load_si512(&a->vec->v[i]);
    g = _mm512_load_si512(&a->vec->v[N/32-1-i]);
    f = _mm512_permutexvar_epi16(permwidx,f);
    g = _mm512_permutexvar_epi16(permwidx,g);
    _mm512_store_si512(&r->vec->v[N/32-1-i],f);
    _mm512_store_si512(&r->vec->v[i],g);
  }
}

void polyvec_sigmam1_ntt(poly *r, const poly *a, size_t len) {
  size_t i;

  for(i=0;i<len;i++)
    poly_sigmam1_ntt(&r[i],&a[i]);
}

void poly_sigma(poly *r, const poly *a, int k) {
  size_t i,j;
  int16_t x;
  poly t;

  j = 0;
  for(i=0;i<N;i++) {
    x = a->vec->c[i];
    x ^= (-(j&N) >> 31) & (x ^ -x);
    t.vec->c[j&(N-1)] = x;
    j += k;
  }

  *r = t;
}

void poly_sigma5(poly *r, const poly *a) {
  poly_sigma(r,a,5);
}

void polyvec_sigma5(poly *r, const poly *a, size_t len) {
  size_t i;

  for(i=0;i<len;i++)
    poly_sigma(&r[i],&a[i],5);
}

void poly_sigma5inv(poly *r, const poly *a) {
  poly_sigma(r,a,3277);
}

void polyvec_sigma5inv(poly *r, const poly *a, size_t len) {
  size_t i;

  for(i=0;i<len;i++)
    poly_sigma(&r[i],&a[i],3277);  // assumes N <= 8192
}

void poly_flip(poly *r, const poly *a) {
  size_t i;
  __m512i f;
  const __m512i ones = _mm512_set1_epi16(1);

  for(i=0;i<N/32;i++) {
    f = _mm512_load_si512(&a->vec->v[i]);
    f = _mm512_sub_epi16(ones,f);
    _mm512_store_si512(&r->vec->v[i],f);
  }
}

void polyvec_flip(poly *r, const poly *a, size_t len) {
  size_t i;

  for(i=0;i<len;i++)
    poly_flip(&r[i],&a[i]);
}

void poly_flip_ntt(poly *r, const poly *a, const pdata *prime) {
  size_t i;
  __m512i f;
  poly t[1];

  f = _mm512_set1_epi16(fpmul(1,prime->s,prime));
  for(i=0;i<N/32;i++)
    _mm512_store_si512(&t->vec->v[i],f);
  poly_ntt(t,t,prime);
  poly_sub(r,t,a);
  poly_sigmam1_ntt(r,r);
}

void polyvec_flip_ntt(poly *r, const poly *a, size_t len, size_t stride, const pdata *prime) {
  size_t i;

  for(i=0;i<len;i++)
    poly_flip_ntt(&r[stride*i],&a[stride*i],prime);
}
