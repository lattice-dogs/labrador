#include <stdint.h>
#include <math.h>
#include "randombytes.h"
#include "gaussian.h"

static const uint64_t rcdt125[11] = {
   4760398266205102531ULL,
   1519614160162702751ULL,
    278740980977854822ULL,
     28213006827475907ULL,
      1542173394262455ULL,
        45012484900334ULL,
          697367627767ULL,
            5716868205ULL,
              24757408ULL,
                 56588ULL,
                    68ULL,
};

/*
 * Compute exp(x) for x such that |x| <= 0.5*ln 2
 *
 * The algorithm used below is derived from the public domain
 * library fdlibm (http://www.netlib.org/fdlibm/e_exp.c).
 *
 */
static double exp_small(double x) {
#define C1 ( 1.66666666666666019037e-01)
#define C2 (-2.77777777770155933842e-03)
#define C3 ( 6.61375632143793436117e-05)
#define C4 (-1.65339022054652515390e-06)
#define C5 ( 4.13813679705723846039e-08)

  double t;
  t = x*x;
  t = x - t*(C1 + t*(C2 + t*(C3 + t*(C4 + t*C5))));
  t = 1 + (x - x*t/(t-2));  // WARNING: Divison not constant time
  return t;

#undef C1
#undef C2
#undef C3
#undef C4
#undef C5
}

static unsigned int rcdtsampler(const uint64_t rcdt[], unsigned int len) {
  unsigned int i,z;
  uint64_t r;

  randombytes((uint8_t *)&r,8);
  r &= (1ULL << 63) - 1;

  z = 0;
  for(i=0;i<len;i++)
    z += (r - rcdt[i]) >> 63;

  return z;
}

static unsigned int BerExp(double x) {
  uint64_t t,u;

  t = (uint64_t)(x*(1/log(2)) + 0.5);
  x -= log(2)*t;
  t = (uint64_t)(exp_small(-x)*exp2(63) + 0.5) >> t;

  randombytes((uint8_t *)&u,8);
  u &= (1ULL << 63) - 1;
  return (u - t) >> 63; // u < t w/ prob t/2^63 = e^-x
}

int64_t gaussian_sampler(double mu, double sigma) {
  int b,k;
  uint8_t bits;
  uint64_t m,mask;
  int64_t r,c1;
  double c0,d,x;

  m = ceil(sigma*(1/1.25));
  mask = m-1;
  mask |= mask >>  1;
  mask |= mask >>  2;
  mask |= mask >>  4;
  mask |= mask >>  8;
  mask |= mask >> 16;
  mask |= mask >> 32;
  do {
    randombytes((uint8_t *)&r,8);
    r &= mask;
  } while((uint64_t)r >= m);

  c0 = (mu+r)*(1.0/m);
  c1 = (int64_t)c0;
  c0 -= c1;  // fractional part
  d = sigma*(1.0/m);
  d = 1/(2*d*d);

  randombytes(&bits,1);
  bits |= 0x80;
  do {
    if(bits <= 1) {
      randombytes(&bits,1);
      bits |= 0x80;
    }
    b = bits & 1;
    bits >>= 1;
    k = rcdtsampler(rcdt125,11);
    k = (-b & (2*k)) - k + b; // bimodal Gaussian
    x = (k-c0)*(k-c0)*d - (k-b)*(k-b)*(1/3.125); // non-negative
  } while(!BerExp(x));

  r = (k+c1)*m - r;
  return r;
}
