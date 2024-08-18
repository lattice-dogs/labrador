#include <stdlib.h>
#include <stdint.h>
#include "codec.h"

/* Encoding of centered Gaussian with standard deviation about 165.74
 * Strategy: As codec2.c but with 8 bit low part
 */

size_t enc(uint8_t *r, int64_t *v, size_t len) {
  size_t i;
  unsigned int off = 0;
  uint8_t *r1 = r+len;
  int8_t l,h;
  uint32_t u;

  r1[0] = 0;
  for(i=0;i<len;i++) {
    r[i] = v[i] & 127;
    l = (r[i] ^ 64) - 64;  // sign extension
    h = (v[i] - l) >> 7;
    l = ~(l >> 7);
    h = (h + l) ^ l;  // flip sign
    l = h >> 7;
    h = (h + l) ^ l;  // absolute value
    h = 2*h + l;
    u = 1 << h;
    r1[0] |= u << off;
    r1[1]  = u >> (8-off);
    r1[2]  = u >> (16-off);
    r1[3]  = u >> (24-off);
    r1 += (off+h+1) >> 3;
    off = (off+h+1) & 7;
  }

  return (r1-r) + (off+7)/8;
}

void dec(int64_t *r, size_t rlen, uint8_t *in, size_t inlen) {
  size_t i;
  unsigned int off = 0;
  uint8_t *in1 = in+rlen;
  int8_t h,l;
  uint8_t b;
  uint32_t u;

  for(i=0;i<rlen;i++) {
    b = in[i] & 127;
    u  = (uint32_t)in1[0] >> off;
    u |= (uint32_t)in1[1] << (8-off);
    u |= (uint32_t)in1[2] << (16-off);
    u |= (uint32_t)in1[3] << (24-off);
    h = __builtin_ctz(u);
    in1 += (off+h+1) >> 3;
    off = (off+h+1) & 7;
    l = -(h & 1);
    h = (h+1) >> 1;
    h = (h + l) ^ l;  // adjust sign
    l = (b ^ 64) - 64;  // sign extension
    r[i] = l;
    l = ~(l >> 7);
    h = (h + l) ^ l;  // flip
    r[i] += (int64_t)h << 7;
  }
}
