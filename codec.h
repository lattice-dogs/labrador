#ifndef CODEC_H
#define CODEC_H

#include <stdlib.h>
#include <stdint.h>

size_t enc(uint8_t *r, int64_t *v, size_t len);
void dec(int64_t *r, size_t rlen, uint8_t *in, size_t inlen);

#endif
