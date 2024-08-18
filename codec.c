#include <stdlib.h>
#include <stdint.h>
#include "codec.h"

/* Encoding of centered Gaussian with standard deviation about 165.74
 * Strategy:
 *     (i) Write h*128 + l for each coefficient with -64 <= l <= 63
 *    (ii) Bitpack all low parts l using 7 bits per coefficients
 *   (iii) Write h = (-1)^sign(h)*abs(h) with sign(h) in {0,1}
 *    (iv) Append sign(h)|1|0...0 (abs(h) zeros) for every coefficient
 */

size_t enc(uint8_t *r, int64_t *v, size_t len) {
  size_t i,j;
  unsigned int off = 0;
  uint8_t *r1 = r+(len*7+7)/8;
  uint8_t b[8],t;
  uint32_t u;

  r1[0] = 0;
  for(i=0;i<len/8;i++) {
    for(j=0;j<8;j++) {
      b[j] = v[8*i+j] & 127;
      t = (v[8*i+j] - (((int64_t)b[j] ^ 64) - 64)) >> 7;  // sign extend
      u = t >> 7;  // sign bit
      t = (t - u) ^ -u;  // absolute value
      u = (u << 1) | 1;
      u = u << t;
      r1[0] |= u << off;
      r1[1]  = u >> (8-off);
      r1[2]  = u >> (16-off);
      r1[3]  = u >> (24-off);
      r1 += (off+t+2) >> 3;
      off = (off+t+2) & 7;
    }

    r[7*i+0]  = b[0];
    r[7*i+0] |= b[1] << 7;
    r[7*i+1]  = b[1] >> 1;
    r[7*i+1] |= b[2] << 6;
    r[7*i+2]  = b[2] >> 2;
    r[7*i+2] |= b[3] << 5;
    r[7*i+3]  = b[3] >> 3;
    r[7*i+3] |= b[4] << 4;
    r[7*i+4]  = b[4] >> 4;
    r[7*i+4] |= b[5] << 3;
    r[7*i+5]  = b[5] >> 5;
    r[7*i+5] |= b[6] << 2;
    r[7*i+6]  = b[6] >> 6;
    r[7*i+6] |= b[7] << 1;
  }

  return (r1-r) + (off+7)/8;
}

void dec(int64_t *r, size_t rlen, uint8_t *in, size_t inlen) {
  size_t i,j;
  unsigned int off = 0;
  uint8_t *in1 = in+(7*rlen+7)/8;
  uint8_t b[8],t;
  uint32_t u;

  for(i=0;i<rlen/8;i++) {
    b[0] = in[0] & 127;
    b[1] = ((in[0] >> 7) | (in[1] << 1)) & 127;
    b[2] = ((in[1] >> 6) | (in[2] << 2)) & 127;
    b[3] = ((in[2] >> 5) | (in[3] << 3)) & 127;
    b[4] = ((in[3] >> 4) | (in[4] << 4)) & 127;
    b[5] = ((in[4] >> 3) | (in[5] << 5)) & 127;
    b[6] = ((in[5] >> 2) | (in[6] << 6)) & 127;
    b[7] = (in[6] >> 1) & 127;
    in += 7;

    for(j=0;j<8;j++) {
      u  = (uint32_t)in1[0] >> off;
      u |= (uint32_t)in1[1] << (8-off);
      u |= (uint32_t)in1[2] << (16-off);
      u |= (uint32_t)in1[3] << (24-off);
      t = __builtin_ctz(u);
      in1 += (off+t+2) >> 3;
      off = (off+t+2) & 7;
      u = (u >> (t+1)) & 1;  // sign bit
      r[8*i+j]  = (((int64_t)t - u) ^ -(int64_t)u) << 7;  // adjust sign
      r[8*i+j] += ((int64_t)b[j] ^ 64) - 64;  // sign extend
    }
  }
}
