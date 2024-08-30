#ifndef DATA_H
#define DATA_H

#define N 64
#ifndef LOGQ
#define LOGQ 32
#endif
#define QBYTES ((LOGQ+7)/8)
#define LOGDELTA log2(1.00444)
#define LIFTS ((128+LOGQ-1)/LOGQ)
#define TAU1 32
#define TAU2 8
#define T 14  // challenge operator norm
#define SLACK 2 // FIXME

#if LOGQ == 24
#define QOFF 3
#define K 5
#define L 2
#define NAMESPACE(s) labrador24_##s
#elif LOGQ == 32
#define QOFF 99
#define K 6
#define L 3
#define NAMESPACE(s) labrador32_##s
#elif LOGQ == 40
#define QOFF 195
#define K 7
#define L 3
#define NAMESPACE(s) labrador40_##s
#elif LOGQ == 48
#define QOFF 59
#define K 8
#define L 4
#define NAMESPACE(s) labrador48_##s
#elif LOGQ == 56
#define QOFF 27
#define K 9
#define L 4
#define NAMESPACE(s) labrador56_##s
#elif LOGQ == 64
#define QOFF 59
#define K 10
#define L 5
#define NAMESPACE(s) labrador64_##s
#else
#error
#endif

#ifndef __ASSEMBLER__

#include <stdint.h>
#include <stddef.h>

#define MIN(a,b) (((a) <= (b)) ? (a) : (b))
#define MAX(a,b) (((a) >= (b)) ? (a) : (b))

typedef struct {
  int16_t p;      // 12 < log2(p) < 14
  int16_t pinv;   // p^-1 mod 2^16
  int16_t v;      // v = round(2^27/p)
  int16_t mont;   // mont = 2^16 mod p
  int16_t montsq; // montsq = 2^32 mod p
  int16_t s;      // s = mont*2^(-14*(L-1)) mod p
  int16_t f;      // f = montsq*2^(14*(L-1)) mod p
  int16_t t;      // t = mont*(P/p)^-1 mod p
  int16_t twist64[64];
  int16_t twist32n[32];
  int16_t twist16[16];
  int16_t twist8n[8];
  int16_t twist4[4];
  int16_t twist2n[2];
  int16_t twist1[2];
  int16_t twist64_pinv[64];
  int16_t twist32n_pinv[32];
  int16_t twist16_pinv[16];
  int16_t twist8n_pinv[8];
  int16_t twist4_pinv[4];
  int16_t twist2n_pinv[2];
  int16_t twist1_pinv[2];
} pdata;

extern const __attribute__((aligned(32))) pdata primes[K];

typedef struct {
  int16_t limbs[L];
} zz;

typedef struct {
  zz q;
  zz pmq;
  zz xvec[K];
} qdata;

static void zz_fromint64(zz *r, int64_t a) {
  size_t i;

  for(i=0;i<L-1;i++) {
    r->limbs[i] = a & 0x3FFF;
    a >>= 14;
  }
  r->limbs[L-1] = a;
}

static int zz_less_than(const zz *a, const zz *b) {
  int i;
  int16_t c = 0;

  for(i=0;i<L;i++)
    c = (a->limbs[i] - b->limbs[i] + c) >> 14;

  return c;
}

static int zz_equal(const zz *a, const zz *b) {
  int i;
  int16_t c = 0;
  int16_t r = 0;

  for(i=0;i<L;i++) {
    c = a->limbs[i] - b->limbs[i] + c;
    r |= c;
    c >>= 14;
  }

  r |= c;
  r &= 0x3FFF;
  r  = -r >> 15;
  r += 1;
  return r;
}

static void zz_add(zz *r, const zz *a, const zz *b) {
  size_t i;
  int16_t c;

  c = 0;
  for(i=0;i<L-1;i++) {
    r->limbs[i] = a->limbs[i] + b->limbs[i] + c;
    c = r->limbs[i] >> 14;
    r->limbs[i] &= 0x3FFF;
  }
  r->limbs[L-1] = a->limbs[L-1] + b->limbs[L-1] + c;
}

static void zz_sub(zz *r, const zz *a, const zz *b) {
  size_t i;
  int16_t c;

  c = 0;
  for(i=0;i<L-1;i++) {
    r->limbs[i] = a->limbs[i] - b->limbs[i] + c;
    c = r->limbs[i] >> 14;
    r->limbs[i] &= 0x3FFF;
  }
  r->limbs[L-1] = a->limbs[L-1] - b->limbs[L-1] + c;
}

__attribute__((aligned(64)))
extern const qdata modulus;

#else
#define _P 0
#define _PINV 1
#define _V 2
#define _MONT 3
#define _MONTSQ 4
#define _S 5
#define _F 6
#define _T 7
#define _TWIST64 8
#define _TWIST32N 72
#define _TWIST16 104
#define _TWIST8N 120
#define _TWIST4 128
#define _TWIST2N 132
#define _TWIST1 134
#define _TWIST64_PINV 136
#define _TWIST32N_PINV 200
#define _TWIST16_PINV 232
#define _TWIST8N_PINV 248
#define _TWIST4_PINV 256
#define _TWIST2N_PINV 260
#define _TWIST1_PINV 262

#endif
#endif
