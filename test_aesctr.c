#include <stdint.h>
#include <stdio.h>
#include <openssl/evp.h>
#include "randombytes.h"
#include "cpucycles.h"
#include "aesctr.h"

static int aes128ctr_evp(uint8_t out[16], const uint8_t key[16], uint64_t nonce, uint64_t ctr)
{
  int i;
  int tmp;
  uint8_t in[16];
  EVP_CIPHER_CTX *ctx;

  ctx = EVP_CIPHER_CTX_new();
  if(!ctx) return -1;

  for(i=0;i<8;i++)
    in[i] = nonce >> 8*i;
  for(i=0;i<8;i++)
    in[8+i] = ctr >> 8*i;

  i = EVP_EncryptInit_ex(ctx, EVP_aes_128_ecb(), 0, key, 0);
  if(i == 1) i = EVP_CIPHER_CTX_set_padding(ctx, 0);
  if(i == 1) i = EVP_EncryptUpdate(ctx, out, &tmp, in, 16);
  if(i == 1) i = EVP_EncryptFinal_ex(ctx, out, &tmp);

  EVP_CIPHER_CTX_free(ctx);
  return (i == 1) ? 0 : -1;
}

int main(void) {
  unsigned int i;
  unsigned long long t[21], overhead;
  __attribute__((aligned(16)))
  uint8_t key[16];
  uint64_t nonce;
  uint64_t ctr = 0;
  size_t nblocks = 10;
  __attribute__((aligned(64)))
  uint8_t out[nblocks*AES128CTR_BLOCKBYTES];
  __attribute__((aligned(64)))
  uint8_t out2[nblocks*AES128CTR_BLOCKBYTES];
  aes128ctr_ctx state;

  overhead = cpucycles_overhead();

  randombytes(key,sizeof(key));
  randombytes((uint8_t*)&nonce,8);

  for(i=0;i<21;i++) {
    t[i] = cpucycles();
    aes128ctr_init(&state,key,nonce);
    aes128ctr_squeezeblocks(out,nblocks,&state);
  }
  for(i=0;i<20;i++)
    printf("aes128ctr: %2d: %llu\n", i, t[i+1] - t[i] - overhead);

  aes128ctr_init(&state,key,nonce);
  aes128ctr_squeezeblocks(out,nblocks,&state);

  for(i=0;i<sizeof(out)/16;i++)
    if(aes128ctr_evp(out2+16*i,key,nonce,ctr++))
      return 1;

  for(i=0;i<sizeof(out);i++)
    if(out2[i] != out[i])
      fprintf(stderr,"ERROR: %d\n", i);

  return 0;
}
