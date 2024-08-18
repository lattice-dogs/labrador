#include <stdint.h>
#include <immintrin.h>
#include <math.h>
#include <stdio.h>
#include "aesctr.h"
#include "data.h"
#include "polx.h"
#include "poly.h"
#include "polz.h"

void polz_print(const polz *a) {
  size_t i,j;

  for(i=0;i<N;i++) {
    printf("%2zu: ",i);
    printf("%5d",a->limbs[0].c[i]);
    for(j=1;j<L;j++)
      printf(" + %5d * 2^%zu",a->limbs[j].c[i],14*j);
    printf("\n");
  }
}

void polzvec_copy(polz *r, const polz *a, size_t len) {
  size_t i,j,k;
  __m512i f;

  for(i=0;i<len;i++) {
    for(j=0;j<L;j++) {
      for(k=0;k<N/32;k++) {
        f = _mm512_load_si512(&a[i].limbs[j].v[k]);
        _mm512_store_si512(&r[i].limbs[j].v[k],f);
      }
    }
  }
}

void polz_getcoeff(zz *r, const polz *a, int k) {
  size_t i;

  for(i=0;i<L;i++)
    r->limbs[i] = a->limbs[i].c[k];
}

void polz_setcoeff(polz *r, const zz *a, int k) {
  size_t i;

  for(i=0;i<L;i++)
    r->limbs[i].c[k] = a->limbs[i];
}

void polz_setcoeff_fromint64(polz *r, int64_t a, int k) {
  size_t i;

  for(i=0;i<L-1;i++) {
    r->limbs[i].c[k] = a & 0x3FFF;
    a >>= 14;
  }
 r->limbs[L-1].c[k] = a;
}

void polzvec_fromint64vec(polz *r, size_t len, size_t deg, const int64_t v[len*deg*N]) {
  size_t i,j,k;

  for(i=0;i<len;i++)
    for(j=0;j<deg;j++)
      for(k=0;k<N;k++)
        polz_setcoeff_fromint64(&r[i*deg+j],v[i*deg*N+k*deg+j],k);
}

int polz_iszero_constcoeff(const polz *a) {
  size_t i;
  int64_t r;

  r = 0;
  for(i=0;i<L;i++)
    r |= (uint16_t)a->limbs[i].c[0];

  r = -r >> 63;
  r += 1;
  return r;
}

// expects centered input
int polz_iszero(const polz *a) {
  size_t i,j;
  int64_t r;

  r = 0;
  for(i=0;i<L;i++)
    for(j=0;j<N;j++)
      r |= (uint16_t)a->limbs[i].c[j];

  r = -r >> 63;
  r += 1;
  return r;
}

int polzvec_iszero(const polz *a, size_t len) {
  size_t i;
  int64_t r;

  r = 0;
  for(i=0;i<len;i++)
    r += polz_iszero(&a[i]);

  r = r-len;
  r >>= 63;
  r += 1;
  return r;
}

static inline size_t uniform(zz *r, const uint8_t *buf) {
  size_t i,j,s,bits;

  i = 0;
  s = 8;
  for(j=0;j<L;j++) {
    bits = (j < L-1) ? 14 : LOGQ - 14*(L-1);
    r->limbs[j] = (int16_t)buf[i] >> (8-s);
    while(s < bits) {
      r->limbs[j] |= (int16_t)buf[++i] << s;
      s += 8;
    }
    s = s - bits;
    r->limbs[j] &= (1 << bits) - 1;
  }

  if(zz_less_than(r,&modulus.q))
    return 1;
  else
    return 0;
}

void polzvec_uniform(polz *r, size_t len, const uint8_t seed[16], uint64_t nonce) {
  size_t i,j,k;
  uint8_t *buf;
  const size_t chunk = (4096/AES128CTR_BLOCKBYTES+QBYTES-1)/QBYTES*AES128CTR_BLOCKBYTES/N;
  size_t nblocks = chunk*N*QBYTES/AES128CTR_BLOCKBYTES;
  aes128ctr_ctx aesctx;
  zz t;

  aes128ctr_init(&aesctx,seed,nonce);

  while(len >= chunk) {
    k = 0;
    buf = (uint8_t*)r + chunk*N*(2*L - QBYTES);
    aes128ctr_squeezeblocks(buf,nblocks,&aesctx);
    for(i=0;i<chunk;i++) {
      for(j=0;j<N;j++) {
        k += uniform(&t,&buf[(i*N+j)*QBYTES]);
        polz_setcoeff(&r[i],&t,j);
      }
    }

    if(k == chunk*N) {
      r += chunk;
      len -= chunk;
    }
  }

  while(len) {
    k = 0;
    nblocks = (len*N*QBYTES+AES128CTR_BLOCKBYTES-1)/AES128CTR_BLOCKBYTES;
    __attribute__((aligned(64)))
    uint8_t buf2[nblocks*AES128CTR_BLOCKBYTES];
    aes128ctr_squeezeblocks(buf2,nblocks,&aesctx);
    for(i=0;i<len;i++) {
      for(j=0;j<N;j++) {
        k += uniform(&t,&buf2[(i*N+j)*QBYTES]);
        polz_setcoeff(&r[i],&t,j);
      }
    }

    if(k == len*N)
      break;
  }
}

// Assumes LOGQ is a multiple of 8
void polz_bitpack(uint8_t r[N*QBYTES], const polz *a) {
  size_t i,j,k,bits;
  __m512i f,g,h;

  k = 0;
  h = _mm512_setzero_si512();
  for(i=0;i<L;i++) {
    bits = (i<L-1) ? 14 : LOGQ-14*(L-1);
    for(j=0;j<N/32;j++) {
      f = _mm512_load_si512(&a->limbs[i].v[j]);
      g = _mm512_slli_epi16(f,k);
      h = _mm512_add_epi16(g,h);
      if(k >= 16-bits) {
        _mm512_storeu_si512(r,h);
        h = _mm512_srli_epi16(f,16-k);
        r += 64;
      }
      k = (k + bits) & 0xF;
    }
  }
}

void polzvec_bitpack(uint8_t *r, const polz *a, size_t len) {
  size_t i;

  for(i=0;i<len;i++)
    polz_bitpack(&r[i*N*QBYTES],&a[i]);
}

void polz_bitunpack(polz *r, const uint8_t buf[(N/32*LOGQ+15)/16*64]) {
  int i,j,k,bits;
  __m512i f,g,h,mask;
  const __m512i mask14 = _mm512_set1_epi16(0x3FFF);
  const __m512i maskhi = _mm512_srli_epi16(mask14,14*L-LOGQ);

  k = 0;
  h = _mm512_setzero_si512();
  for(i=0;i<L;i++) {
    if(i < L-1) {
      bits = 14;
      mask = mask14;
    }
    else {
      bits = LOGQ - 14*(L-1);
      mask = maskhi;
    }

    for(j=0;j<N/32;j++) {
      if(k < bits) {
        f = _mm512_load_si512((__m512i*)buf);
        g = _mm512_slli_epi16(f,k);
        h = _mm512_add_epi16(g,h);
        g = _mm512_and_si512(h,mask);
        h = _mm512_srli_epi16(f,bits-k);
        buf += 64;
      }
      else {
        g = _mm512_and_si512(h,mask);
        h = _mm512_srli_epi16(h,bits);
      }
      _mm512_store_si512(&r->limbs[i].v[j],g);
      k = (k - bits) & 0xF;
    }
  }
}

void polzvec_almostuniform(polz *r, size_t len, const uint8_t seed[16], uint64_t nonce) {
  size_t i;
  uint8_t *buf;
  size_t chunk,nblocks;
  aes128ctr_ctx aesctx;

#if LOGQ%4 == 0
  chunk = (4*4096/AES128CTR_BLOCKBYTES+LOGQ-1)/LOGQ*2*AES128CTR_BLOCKBYTES/N;
#elif
  chunk = (4096/AES128CTR_BLOCKBYTES+LOGQ-1)/LOGQ*8*AES128CTR_BLOCKBYTES/N;
#endif

  nblocks = chunk*N/AES128CTR_BLOCKBYTES*LOGQ/8;
  aes128ctr_init(&aesctx,seed,nonce);

  while(len >= chunk) {
    buf = (uint8_t*)r + chunk*(2*N*L - N*LOGQ/8);
    aes128ctr_squeezeblocks(buf,nblocks,&aesctx);
    for(i=0;i<chunk;i++)
      polz_bitunpack(&r[i],&buf[i*N*LOGQ/8]);
    len -= chunk;
    r += chunk;
  }

  if(len) {
    nblocks = (len*N*LOGQ/8+AES128CTR_BLOCKBYTES-1)/AES128CTR_BLOCKBYTES;
    __attribute__((aligned(64)))
    uint8_t buf2[nblocks*AES128CTR_BLOCKBYTES];
    aes128ctr_squeezeblocks(buf2,nblocks,&aesctx);
    for(i=0;i<len;i++)
      polz_bitunpack(&r[i],&buf2[i*N*LOGQ/8]);
  }
}

double polzvec_norm(const polz *a, size_t len) {
  size_t i,j,k;
  long double r,t;
  long double q,hq;

  q = ldexpl(1,LOGQ) - QOFF;
  hq = (q-1)/2;

  r = 0;
  for(i=0;i<len;i++) {
    for(j=0;j<N;j++) {
      t = 0;
      for(k=0;k<L;k++)
        t += ldexpl(a[i].limbs[k].c[j],14*k);
      if(t > hq) t -= q;
      r += t*t;
    }
  }

  return sqrtl(r);
}

void polz_reduce(polz *r) {
  size_t i,j;
  __m512i f,g,h;
  const __m512i qoff = _mm512_set1_epi16(QOFF);
  const __m512i qofft4 = _mm512_slli_epi16(qoff,2);
  const __m512i mask14 = _mm512_set1_epi16(0x3FFF);
  const __m512i maskhi = _mm512_srli_epi16(mask14,14*L-LOGQ);

  for(i=0;i<N/32;i++) {
    f = _mm512_load_si512(&r->limbs[L-1].v[i]);
    g = _mm512_srai_epi16(f,LOGQ-14*(L-1));
    f = _mm512_and_si512(f,maskhi);
    _mm512_store_si512(&r->limbs[L-1].v[i],f);

    f = _mm512_load_si512(&r->limbs[0].v[i]);
    h = _mm512_mullo_epi16(g,qoff);
    g = _mm512_mulhi_epi16(g,qofft4);
    h = _mm512_and_si512(h,mask14);
    f = _mm512_add_epi16(f,h);
    h = _mm512_srai_epi16(f,14);
    f = _mm512_and_si512(f,mask14);
    g = _mm512_add_epi16(g,h);
    _mm512_store_si512(&r->limbs[0].v[i],f);

    for(j=1;j<L;j++) {
      f = _mm512_load_si512(&r->limbs[j].v[i]);
      f = _mm512_add_epi16(f,g);
      if(j<L-1) {
        g = _mm512_srai_epi16(f,14);
        f = _mm512_and_si512(f,mask14);
      }
      _mm512_store_si512(&r->limbs[j].v[i],f);
    }
  }
}

void polzvec_reduce(polz *r, size_t len) {
  size_t i;

  for(i=0;i<len;i++)
    polz_reduce(&r[i]);
}

void polz_caddq(polz *r) {
  size_t i,j;
  __m512i f,c;
  __mmask32 mask;
  __m512i qq[L];
  const __m512i mask14 = _mm512_set1_epi16(0x3FFF);

  for(i=0;i<L;i++)
    qq[i] = _mm512_set1_epi16(modulus.q.limbs[i]);

  for(i=0;i<N/32;i++) {
    f = _mm512_load_si512(&r->limbs[L-1].v[i]);
    mask = _mm512_movepi16_mask(f);
    for(j=0;j<L;j++) {
      f = _mm512_load_si512(&r->limbs[j].v[i]);
      f = _mm512_mask_add_epi16(f,mask,f,qq[j]);
      if(j > 0) f = _mm512_add_epi16(f,c);
      if(j < L-1) {
        c = _mm512_srai_epi16(f,14);
        f = _mm512_and_si512(f,mask14);
      }
      _mm512_store_si512(&r->limbs[j].v[i],f);
    }
  }
}

void polzvec_caddq(polz *r, size_t len) {
  size_t i;

  for(i=0;i<len;i++)
    polz_caddq(&r[i]);
}

void polz_center(polz *r) {
  size_t i,j;
  __m512i f,c;
  __mmask32 mask;
  __m512i qq[L], hq[L];
  const __m512i mask14 = _mm512_set1_epi16(0x3FFF);

  for(i=0;i<L;i++)
    qq[i] = _mm512_set1_epi16(modulus.q.limbs[i]);
  for(i=0;i<L-1;i++) {
    hq[i] = _mm512_srli_epi16(qq[i],1);
    f = _mm512_slli_epi16(qq[i+1],13);
    hq[i] = _mm512_add_epi16(hq[i],f);
    hq[i] = _mm512_and_si512(hq[i],mask14);
  }
  hq[L-1] = _mm512_srli_epi16(qq[L-1],1);

  for(i=0;i<N/32;i++) {
    for(j=0;j<L;j++) {
      f = _mm512_load_si512(&r->limbs[j].v[i]);
      f = _mm512_sub_epi16(hq[j],f);
      if(j > 0) f = _mm512_add_epi16(f,c);
      if(j < L-1) c = _mm512_srai_epi16(f,14);
    }
    mask = _mm512_movepi16_mask(f);

    for(j=0;j<L;j++) {
      f = _mm512_load_si512(&r->limbs[j].v[i]);
      f = _mm512_mask_sub_epi16(f,mask,f,qq[j]);
      if(j > 0) f = _mm512_add_epi16(f,c);
      if(j < L-1) {
        c = _mm512_srai_epi16(f,14);
        f = _mm512_and_si512(f,mask14);
      }
      _mm512_store_si512(&r->limbs[j].v[i],f);
    }
  }
}

void polzvec_center(polz *r, size_t len) {
  size_t i;

  for(i=0;i<len;i++)
    polz_center(&r[i]);
}

static void polz_topoly_montgomery(poly *r, const polz *a, const pdata *prime) {
  size_t i,j;
  __m512i f,g,h;
  const __m512i pinvt4 = _mm512_slli_epi16(_mm512_set1_epi16(prime->pinv),2);
  const __m512i p = _mm512_set1_epi16(prime->p);

  for(i=0;i<N/32;i++) {
    g = _mm512_load_si512(&a->limbs[0].v[i]);
    for(j=1;j<L;j++) {
      h = _mm512_mullo_epi16(g,pinvt4);
      f = _mm512_load_si512(&a->limbs[j].v[i]);
      g = _mm512_srai_epi16(g,14);
      h = _mm512_mulhi_epu16(h,p);
      g = _mm512_add_epi16(f,g);
      g = _mm512_sub_epi16(g,h);
    }
    _mm512_store_si512(&r->vec->v[i],g);
  }

  //Scaling factor multiplied during NTT
  //poly_scale(r,r,(1LL << (14*(L-1)+16)) % prime->p,prime);
}

static void polzvec_topolyvec_montgomery(poly *r, const polz *a, size_t len, size_t stride, const pdata *prime) {
  size_t i;

  for(i=0;i<len;i++)
    polz_topoly_montgomery(&r[stride*i],&a[i],prime);
}

void polz_topoly(poly *r, const polz *a) {
  size_t i;
  __m512i f;

  for(i=0;i<N/32;i++) {
    f = _mm512_load_si512(&a->limbs[0].v[i]);
    _mm512_store_si512(&r->vec->v[i],f);
  }
}

void polzvec_topolyvec(poly *r, const polz *a, size_t len, size_t stride) {
  size_t i;

  for(i=0;i<len;i++)
    polz_topoly(&r[stride*i],&a[i]);
}

void polz_frompoly(polz *r, const poly *a) {
  size_t i,j;
  __m512i f;
  const __m512i zero = _mm512_setzero_si512();

  for(i=0;i<N/32;i++) {
    f = _mm512_load_si512(&a->vec->v[i]);
    _mm512_store_si512(&r->limbs[0].v[i],f);
  }

  for(i=1;i<L;i++)
    for(j=0;j<N/32;j++)
      _mm512_store_si512(&r->limbs[i].v[j],zero);
}

void polzvec_frompolyvec(polz *r, const poly *a, size_t len) {
  size_t i;

  for(i=0;i<len;i++)
    polz_frompoly(&r[i],&a[i]);
}

void polz_topolx(polx *r, const polz *a) {
  size_t i;

  for(i=0;i<K;i++)
    polz_topoly_montgomery(&r->vec[i],a,&primes[i]);

  polx_ntt(r,r);
}

void polzvec_topolxvec(polx *r, const polz *a, size_t len) {
  size_t i;

  for(i=0;i<K;i++)
    polzvec_topolyvec_montgomery(&r->vec[i],a,len,K,&primes[i]);

  polxvec_ntt(r,r,len);
}

static void zz_poly_mul(polz *r, const zz *a, const poly *b) {
  size_t i,j;
  int bits;
  __m512i e[L],f,ff;
  __m512i g,h,k,l;
  __m512i mask,qoff;

  for(i=0;i<L;i++)
    e[i] =  _mm512_set1_epi16(a->limbs[i]);  // 2^bits

  for(i=0;i<N/32;i++) {
    f = _mm512_load_si512(&b->vec->v[i]);  // 2^14
    ff = _mm512_slli_epi16(f,2);  // 2^16

    bits = LOGQ - 14*(L-1);
    mask = _mm512_set1_epi16((1 << bits) - 1);
    h = _mm512_mullo_epi16(f,e[L-1]);  // 2^16
    g = _mm512_mulhi_epu16(f,e[L-1]);  // 2^(bits-2)
    k = _mm512_srli_epi16(h,bits);  // 2^(16-bits)
    g = _mm512_slli_epi16(g,16-bits);  // 2^14
    h = _mm512_and_si512(h,mask);  // 2^bits
    k = _mm512_add_epi16(k,g);  // 2^14
    _mm512_store_si512(&r->limbs[L-1].v[i],h);  // 2^bits

    // Assumes L > 1
    bits = 14;
    mask = _mm512_set1_epi16(0x3FFF);
    qoff = _mm512_set1_epi16(QOFF);  // 2^13
    h = _mm512_mullo_epi16(f,e[0]);  // 2^16
    l = _mm512_mullo_epi16(qoff,k);  // 2^16
    qoff = _mm512_slli_epi16(qoff,2);  // 2^15
    g = _mm512_mulhi_epu16(ff,e[0]);  // 2^14
    k = _mm512_mulhi_epu16(qoff,k);  // 2^13
    h = _mm512_and_si512(h,mask);  // 2^14
    l = _mm512_and_si512(l,mask);  // 2^14
    h = _mm512_add_epi16(h,l);  // 2^15-1
    l = _mm512_srli_epi16(h,bits);  // 2
    h = _mm512_and_si512(h,mask);  // 2^14
    k = _mm512_add_epi16(k,g);  // 2^14+2^13-1
    k = _mm512_add_epi16(k,l);  // 2^14+2^13
    _mm512_store_si512(&r->limbs[0].v[i],h);  // 2^14

    for(j=1;j<L-1;j++) {
      h = _mm512_mullo_epi16(f,e[j]);  // 2^16
      g = _mm512_mulhi_epu16(ff,e[j]);  // 2^13
      h = _mm512_and_si512(h,mask);  // 2^14
      h = _mm512_add_epi16(h,k);  // 2^15+2^13-2
      k = _mm512_srli_epi16(h,bits);  // 2+1
      h = _mm512_and_si512(h,mask);  // 2^14
      k = _mm512_add_epi16(k,g);  // 2^13+2
      _mm512_store_si512(&r->limbs[j].v[i],h);  // 2^14
    }

    g = _mm512_load_si512(&r->limbs[L-1].v[i]);  // 2^bits
    g = _mm512_add_epi16(g,k);  // 2^bits+2^14+2^13-1
    _mm512_store_si512(&r->limbs[L-1].v[i],g);  // 2^bits+2^14+2^13-1
  }
}

static void zz_poly_fma(polz *r, const zz *a, const poly *b) {
  size_t i,j;
  int bits;
  __m512i e[L],f,ff;
  __m512i g,h,k,l;
  __m512i mask,qoff;

  for(i=0;i<L;i++)
    e[i] = _mm512_set1_epi16(a->limbs[i]);  // 2^bits

  for(i=0;i<N/32;i++) {
    f = _mm512_load_si512(&b->vec->v[i]);  // 2^14 (p_i)
    ff = _mm512_slli_epi16(f,2);  // 2^16

    bits = LOGQ - 14*(L-1);
    mask = _mm512_set1_epi16((1 << bits) - 1);
    h = _mm512_mullo_epi16(f,e[L-1]);  // 2^16
    g = _mm512_mulhi_epu16(f,e[L-1]);  // 2^(bits-2)
    k = _mm512_srli_epi16(h,bits);  // 2^(16-bits)
    g = _mm512_slli_epi16(g,16-bits);  // 2^14
    h = _mm512_and_si512(h,mask);  // 2^bits
    g = _mm512_add_epi16(g,k);  // 2^14
    k = _mm512_load_si512(&r->limbs[L-1].v[i]);  // 2^bits+2^15
    h = _mm512_add_epi16(h,k);  // 2^(bits+1)+2^15-1
    k = _mm512_srli_epi16(h,bits);  // 2^(15-bits)+2
    h = _mm512_and_si512(h,mask);  // 2^bits
    k = _mm512_add_epi16(k,g);  // 2^14+2^12+1
    _mm512_store_si512(&r->limbs[L-1].v[i],h);  // 2^bits

    bits = 14;
    qoff = _mm512_set1_epi16(QOFF);  // 2^13
    mask = _mm512_set1_epi16(0x3FFF);

    h = _mm512_mullo_epi16(f,e[0]);  // 2^16
    l = _mm512_mullo_epi16(qoff,k);  // 2^16
    qoff = _mm512_slli_epi16(qoff,2);  // 2^15
    g = _mm512_mulhi_epu16(ff,e[0]);  // 2^14
    k = _mm512_mulhi_epu16(qoff,k);  // 2^14
    h = _mm512_and_si512(h,mask);  // 2^14
    l = _mm512_and_si512(l,mask);  // 2^14
    h = _mm512_add_epi16(h,l);  // 2^15-1
    l = _mm512_load_si512(&r->limbs[0].v[i]);  // 2^14
    h = _mm512_add_epi16(h,l);  // 2^15+2^14-2
    l = _mm512_srli_epi16(h,bits);  // 2+1
    h = _mm512_and_si512(h,mask);  // 2^14
    k = _mm512_add_epi16(k,g);  // 2^15-1
    k = _mm512_add_epi16(k,l);  // 2^15+1
    _mm512_store_si512(&r->limbs[0].v[i],h);  // 2^14

    for(j=1;j<L-1;j++) {
      h = _mm512_mullo_epi16(f,e[j]);  // 2^16
      g = _mm512_mulhi_epu16(ff,e[j]);  // 2^14
      h = _mm512_and_si512(h,mask);  // 2^14
      h = _mm512_add_epi16(h,k);  // 2^15+2^14
      k = _mm512_load_si512(&r->limbs[j].v[i]);  // 2^14
      h = _mm512_add_epi16(h,k);  // 2^16-1
      k = _mm512_srli_epi16(h,bits);  // 2^2
      h = _mm512_and_si512(h,mask);  // 2^14
      k = _mm512_add_epi16(k,g);  // 2^14+2^2-1
      _mm512_store_si512(&r->limbs[j].v[i],h);  // 2^14
    }

    g = _mm512_load_si512(&r->limbs[L-1].v[i]);  // 2^bits
    g = _mm512_add_epi16(g,k);  // 2^bits+2^15
    _mm512_store_si512(&r->limbs[L-1].v[i],g);  // 2^bits+2^15
  }
}

/* Explicit CRT mod q:
 * a_i = a mod p_i, |a| <= (P-1)/2
 * t_i = (P/p_i)^-1 mod p_i
 * alpha_i = a_it_i mod p_i
 * Explicit CRT: a = (\sum_i alpha_i/p_i - round(\sum_i alpha_i/p_i))P
 * Explicit CRT mod q: \sum_i alpha_i(P/p_i mod q) - round(\sum_i alpha_i/p_i)(P mod q)
 */
void polz_frompolx(polz *r, const polx *a) {
  size_t i;
  polx b;
  poly k = {0};
  __m512i f;
  const __m512i shift = _mm512_set1_epi16(1024);

  polx_invntt(&b,a);

  for(i=0;i<K;i++) {
    poly_scale(&b.vec[i],&b.vec[i],primes[i].t,&primes[i]);  // alpha_i
    poly_caddp(&b.vec[i],&primes[i]);  // needed for higher precision
    poly_quot_add(&k,&b.vec[i],&primes[i]);  // TODO: Overflow possible for K > 8
  }

  for(i=0;i<N/32;i++) {
    f = _mm512_load_si512(&k.vec->v[i]);
    f = _mm512_add_epi16(f,shift);
    f = _mm512_srli_epi16(f,11);  // round(\sum_i alpha_i/p_i)
    _mm512_store_si512(&k.vec->v[i],f);
  }

  zz_poly_mul(r,&modulus.pmq,&k);
  for(i=0;i<K;i++)
    zz_poly_fma(r,&modulus.xvec[i],&b.vec[i]);

  polz_reduce(r);
}

static void frompolxvec(polz *r, const polx *a, size_t len) {
  size_t i,j;
  polx b[len];
  poly k[len];
  __m512i f;
  const __m512i shift = _mm512_set1_epi16(1024);

  polxvec_invntt(b,a,len);

  polyvec_setzero(k,len);
  for(i=0;i<K;i++) {
    polyvec_scale(&b->vec[i],&b->vec[i],len,K,primes[i].t,&primes[i]);  // alpha_i
    polyvec_caddp(&b->vec[i],len,K,&primes[i]);  // needed for higher precision
    polyvec_quot_add(k,&b->vec[i],len,K,&primes[i]);  // TODO: Overflow possible for K > 8
  }

  for(i=0;i<len;i++) {
    for(j=0;j<N/32;j++) {
      f = _mm512_load_si512(&k[i].vec->v[j]);
      f = _mm512_add_epi16(f,shift);
      f = _mm512_srli_epi16(f,11);  // round(\sum_i alpha_i/p_i)
      _mm512_store_si512(&k[i].vec->v[j],f);
    }
  }

  for(i=0;i<len;i++) {
    zz_poly_mul(&r[i],&modulus.pmq,&k[i]);
    for(j=0;j<K;j++)
      zz_poly_fma(&r[i],&modulus.xvec[j],&b[i].vec[j]);
  }

  polzvec_reduce(r,len);
}

void polzvec_frompolxvec(polz *r, const polx *a, size_t len) {
  while(len >= 16) {
    frompolxvec(r,a,16);
    r += 16;
    a += 16;
    len -= 16;
  }

  frompolxvec(r,a,len);
}

void polz_add(polz *r, const polz *a, const polz *b) {
  size_t i,j;
  __m512i f,g,c;
  const __m512i mask14 = _mm512_set1_epi16(0x3FFF);

  for(i=0;i<N/32;i++) {
    for(j=0;j<L;j++) {
      f = _mm512_load_si512(&a->limbs[j].v[i]);
      g = _mm512_load_si512(&b->limbs[j].v[i]);
      f = _mm512_add_epi16(f,g);
      if(j > 0) f = _mm512_add_epi16(f,c);
      if(j < L-1) {
        c = _mm512_srai_epi16(f,14);
        f = _mm512_and_si512(f,mask14);
      }
      _mm512_store_si512(&r->limbs[j].v[i],f);
    }
  }
}

void polzvec_add(polz *r, const polz *a, const polz *b, size_t len) {
  size_t i;

  for(i=0;i<len;i++)
    polz_add(&r[i],&a[i],&b[i]);
}

void polz_sub(polz *r, const polz *a, const polz *b) {
  size_t i,j;
  __m512i f,g,c;
  const __m512i mask14 = _mm512_set1_epi16(0x3FFF);

  for(i=0;i<N/32;i++) {
    for(j=0;j<L;j++) {
      f = _mm512_load_si512(&a->limbs[j].v[i]);
      g = _mm512_load_si512(&b->limbs[j].v[i]);
      f = _mm512_sub_epi16(f,g);
      if(j > 0) f = _mm512_add_epi16(f,c);
      if(j < L-1) {
        c = _mm512_srai_epi16(f,14);
        f = _mm512_and_si512(f,mask14);
      }
      _mm512_store_si512(&r->limbs[j].v[i],f);
    }
  }
}

void polzvec_sub(polz *r, const polz *a, const polz *b, size_t len) {
  size_t i;

  for(i=0;i<len;i++)
    polz_sub(&r[i],&a[i],&b[i]);
}

void polz_slli(polz *r, const polz *a, int s) {
  size_t i,j;
  __m512i f,g,h;
  const __m512i mask14 = _mm512_set1_epi16(0x3FFF);

  for(i=0;i<N/32;i++) {
    for(j=0;j<L;j++) {
      f = _mm512_load_si512(&a->limbs[j].v[i]);
      g = _mm512_slli_epi16(f,s);
      if(j) g = _mm512_add_epi16(g,h);
      if(j<L-1) {
        g = _mm512_and_si512(g,mask14);
        h = _mm512_srai_epi16(f,14-s);
      }
      _mm512_store_si512(&r->limbs[j].v[i],g);
    }
  }
}

void polzvec_slli(polz *r, const polz *a, size_t len, int s) {
  size_t i;

  for(i=0;i<len;i++)
    polz_slli(&r[i],&a[i],s);
}

void polz_mul(polz *r, const polz *a, const polz *b) {
  polx f,g;

  polz_topolx(&f,a);
  polz_topolx(&g,b);
  polx_mul(&f,&f,&g);
  polz_frompolx(r,&f);
}

void polz_poly_mul(polz *r, const polz *a, const poly *b) {
  polx f,g;

  polz_topolx(&f,a);
  polx_frompoly(&g,b);
  polx_mul(&f,&f,&g);
  polz_frompolx(r,&f);
}

// expects centered input
void polz_split(poly *lo, polz *hi, const polz *a, size_t d) {
  size_t i,j;
  __m512i f,g,h;
  const __m512i mask14 = _mm512_set1_epi16(0x3FFF);

  for(i=0;i<N/32;i++) {
    f = _mm512_load_si512(&a->limbs[0].v[i]);
    g = _mm512_slli_epi16(f,16-d);
    g = _mm512_srai_epi16(g,16-d);  // mod 2^d sign extended
    _mm512_store_si512(&lo->vec->v[i],g);
    f = _mm512_sub_epi16(f,g);  // zero mod 2^d
    f = _mm512_srai_epi16(f,d);
    for(j=0;j<L-1;j++) {
      g = _mm512_load_si512(&a->limbs[j+1].v[i]);
      h = _mm512_slli_epi16(g,14-d);
      h = _mm512_and_si512(h,mask14);
      f = _mm512_add_epi16(f,h);
      h = _mm512_srai_epi16(f,14);  // carry
      f = _mm512_and_si512(f,mask14);
      _mm512_store_si512(&hi->limbs[j].v[i],f);
      f = _mm512_srai_epi16(g,d);
      f = _mm512_add_epi16(f,h);
    }
    _mm512_store_si512(&hi->limbs[L-1].v[i],f);
  }
}

void polzvec_split(poly *lo, polz *hi, const polz *a, size_t len, size_t d) {
  size_t i;

  for(i=0;i<len;i++)
    polz_split(&lo[i],&hi[i],&a[i],d);
}

// expects centered input
void polz_decompose(poly *r, const polz *a, size_t stride, size_t t, size_t d) {
  size_t i,j,k,s;
  __m512i f,g,h;
  const __m512i mask = _mm512_set1_epi16((1<<d) - 1);

  for(i=0;i<N/32;i++) {
    f = _mm512_load_si512(&a->limbs[0].v[i]);
    k = 1;
    s = 14;
    for(j=0;j<t-1;j++) {
      if(s < d) {
        g = _mm512_load_si512(&a->limbs[k++].v[i]);
        h = _mm512_slli_epi16(g,s);
        h = _mm512_and_si512(h,mask);
        f = _mm512_add_epi16(f,h);
        h = _mm512_slli_epi16(f,16-d);
        h = _mm512_srai_epi16(h,16-d);  // mod 2^d sign extended
        _mm512_store_si512(&r[stride*j].vec->v[i],h);
        f = _mm512_sub_epi16(f,h);  // zero mod 2^d
        f = _mm512_srai_epi16(f,d);
        g = _mm512_srai_epi16(g,d-s);
        f = _mm512_add_epi16(f,g);
        s += 14-d;
      }
      else {
        g = _mm512_slli_epi16(f,16-d);
        g = _mm512_srai_epi16(g,16-d);  // mod 2^d sign extended
        _mm512_store_si512(&r[stride*j].vec->v[i],g);
        f = _mm512_sub_epi16(f,g);  // zero mod 2^d
        f = _mm512_srai_epi16(f,d);
        s -= d;
      }
    }

    if(k < L) {
      g = _mm512_load_si512(&a->limbs[k++].v[i]);
      g = _mm512_slli_epi16(g,s);
      f = _mm512_add_epi16(f,g);
    }
    _mm512_store_si512(&r[stride*(t-1)].vec->v[i],f);
  }
}

void polzvec_decompose(poly *r, const polz *a, size_t len, size_t t, size_t d) {
  size_t i;

  if(t==1) {
    polzvec_topolyvec(r,a,len,1);
    return;
  }

  for(i=0;i<len;i++)
    polz_decompose(&r[i],&a[i],len,t,d);
}

void polz_decompose_topolx(polx *r, const polz *a, size_t stride, size_t t, size_t d) {
  size_t i;
  poly b[t];

  polz_decompose(b,a,1,t,d);
  for(i=0;i<t;i++)
    polx_frompoly(&r[stride*i],&b[i]);
}

void polzvec_decompose_topolxvec(polx *r, const polz *a, size_t len, size_t stride, size_t t, size_t d) {
  size_t i;
  poly b[16*t];

  if(len < stride)
    for(i=0;i<t;i++)
      polxvec_setzero(&r[stride*i+len],stride-len);

  while(len >= 16) {
    polzvec_decompose(b,a,16,t,d);
    for(i=0;i<t;i++)
      polxvec_frompolyvec(&r[stride*i],&b[16*i],16);
    r += 16;
    a += 16;
    len -= 16;
  }

  if(len) {
    polzvec_decompose(b,a,len,t,d);
    for(i=0;i<t;i++)
      polxvec_frompolyvec(&r[stride*i],&b[len*i],len);
  }
}

void polz_reconstruct(polz *r, const poly *a, size_t stride, size_t t, size_t d) {
  size_t i,j,k,s;
  __m512i f,g,h;
  const __m512i mask = _mm512_set1_epi16(0x3FFF);

  for(i=0;i<N/32;i++) {
    h = _mm512_setzero_si512();
    k = s = 0;
    for(j=0;j<L;j++) {
      while(k < t && (j == L-1 || s < 14-d)) {
        f = _mm512_load_si512(&a[stride*k++].vec->v[i]);
        f = _mm512_slli_epi16(f,s);
        h = _mm512_add_epi16(h,f);
        s += d;
      }  // k == t || (j < L-1 && 14-d <= s < 14)
      if(k < t) {
        f = _mm512_load_si512(&a[stride*k++].vec->v[i]);
        g = _mm512_slli_epi16(f,s);
        g = _mm512_and_si512(g,mask);
        h = _mm512_add_epi16(h,g);
        g = _mm512_and_si512(h,mask);
        _mm512_store_si512(&r->limbs[j].v[i],g);
        h = _mm512_srai_epi16(h,14);
        g = _mm512_srai_epi16(f,14-s);
        h = _mm512_add_epi16(h,g);
        s -= 14-d;
      }  // k == t || 0 <= s < d
      else if(j < L-1) {
        g = _mm512_and_si512(h,mask);
        _mm512_store_si512(&r->limbs[j].v[i],g);
        h = _mm512_srai_epi16(h,14);
      }
      else
        _mm512_store_si512(&r->limbs[j].v[i],h);
    }
  }
}

void polzvec_reconstruct(polz *r, const poly *a, size_t len, size_t t, size_t d) {
  size_t i;

  if(t==1) {
    polzvec_frompolyvec(r,a,len);
    return;
  }

  for(i=0;i<len;i++)
    polz_reconstruct(&r[i],&a[i],len,t,d);
}

void polz_sigmam1(polz *r, const polz *a) {
  polyvec_sigmam1((poly*)r,(poly*)a,L);
}

void polzvec_sigmam1(polz *r, const polz *a, size_t len) {
  polyvec_sigmam1((poly*)r,(poly*)a,L*len);
}
