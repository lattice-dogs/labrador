#!/usr/bin/env python3

import sys
import math

if len(sys.argv) == 2:
  logq = int(sys.argv[1])
else:
  logq = 32

primes = [
  [ 7937, 29],
  [ 7681,202],
  [ 7297,131],
  [ 6529,118],
  [ 4993, 66],
  [ 4481, 42],
#  [16001,523],
#  [15361,100],
#  [15233, 42],
#  [14593,491],
#  [14081,657],
  [13697,175],
  [13441,110],
  [13313,272],
  [12289, 81],
  [12161,182],
  [11777,309],
  [11393,182],
  [10753, 94],
  [10369,171],
  [ 9857, 58],
  [ 9601, 31],
  [ 9473, 28]
]

qoffs = [
     -1,   -1,   -1,    3,    3,    3,    3,   19,
     27,    3,    3,   19,    3,   75,    3,   19,
     99,   91,   11,   19,    3,   19,    3,   27,
      3,   91,   27,  115,  299,    3,   35,   19,
     99,  355,  131,  451,  243,  123,  107,   19,
    195,   75,   11,   67,  539,  139,  635,  115,
     59,  123,   27,  139,  395,  315,  131,   67,
     27,  195,   27,   99,  107,  259,  171,  259,
     59,  115,  203,   19,   83,   19,   35,  411,
    107,  475,   35,  427,  123,   43,   11,   67,
   1307,   51,  315,  139,   35,   19,   35,   67,
    299,   99,   75,  315,   83,   51,    3,  211,
    147,  595,   51,  115,   99,   99,  483,  339,
    395,  139, 1187,  171,   59,   91,  195,  835,
     75,  211,   11,   67,    3,  451,  563,  867,
    395,  531,    3,   67,   59,  579,  203,  507,
    275,  315,   27,  315,  347,   99,  603,  795,
    243,  339,  203,  187,   27,  171, 1491,  355,
     83,  355, 1371,  387,  347,   99,    3,  195,
    539,  171,  243,  499,  195,   19,  155,   91,
     75, 1011,  627,  867,  155,  115, 1811,  771,
   1467,  643,  195,   19,  155,  531,    3,  267,
    563,  339,  563,  507,  107,  283,  267,  147,
     59,  339,  371, 1411,  363,  819,   11,   19,
    915,  123,   75,  915,  459,   75,  627,  459,
     75, 1035,  195,  187, 1515, 1219, 1443,   91,
    299,  451,  171, 1099,   99,    3,  395, 1147,
    683,  675,  243,  355,  395,    3,  875,  235,
    363, 1131,  155,  835,  723,   91,   27,  235,
    875,    3,   83,  259,  875, 1515,  731,  531,
    467,  819,  267,  475, 1923,  163,  107,  411,
    387,   75, 2331,  355, 1515, 1723, 1427,   19
]

def centermod(a,b):
  r = a % abs(b)
  if r > (b-1)//2: r -= b
  return r

def bitrev6(a):
  t  = (a &  1) << 5
  t |= (a &  2) << 3
  t |= (a &  4) << 1
  t |= (a &  8) >> 1
  t |= (a & 16) >> 3
  t |= (a & 32) >> 5
  return t

q = 2**logq - qoffs[logq]
nlimbs = math.ceil(logq/14)
P = 1
nprimes = 0
for prime in primes:
  P *= prime[0]
  nprimes += 1
  if P > 128*q**2: break # FIXME: Map to [0,P-1] in CRT

print("#include <stdint.h>")
print("#include \"data.h\"")
print()

print("#define N 64")
print("#define LOGQ %d"%logq)
print("#define QOFF %d"%qoffs[logq])
print("#define K %d"%nprimes)
print("#define L %d"%nlimbs)
print()

print("__attribute__((aligned(64)))")
print("const pdata primes[%d] = {"%nprimes)

for prime in primes[:nprimes]:
  p = prime[0]
  zeta = prime[1]
  pinv = centermod(pow(p,-1,2**16),2**16)
  v = math.floor(2**27/p + 0.5)
  mont = centermod(2**16,p)
  montsq = centermod(2**32,p)
#  s = centermod(mont*pow(2**(14*(nlimbs-1)),-1,p),p)
#  f = centermod(mont*2**5,p) # 2*5-16+6=0
#  t = centermod(mont*pow(P//p,-1,p)*2**(2*14*(nlimbs-1)),p)
  s = centermod(mont*pow(2,-14*(nlimbs-1),p),p)
  f = centermod(montsq*pow(2,14*(nlimbs-1),p),p)
  t = centermod(pow(64*P//p,-1,p),p)

  print("  {")
  print("    .p = %d,"%p)
  print("    .pinv = %d,"%pinv)
  print("    .v = %d,"%v)
  print("    .mont = %d,"%mont)
  print("    .montsq = %d,"%montsq)
  print("    .s = %d,"%s)
  print("    .f = %d,"%f)
  print("    .t = %d,"%t)

#  print("    .zetas = {")
#  for i in range(64//8):
#    print("      ",end='')
#    for j in range(8):
#      print("%6d,"%centermod(mont*pow(zeta,bitrev6(8*i+j),p),p),end='')
#    print()
#  print("    },")
#  print("    .zetas_pinv = {")
#  for i in range(64//8):
#    print("      ",end='')
#    for j in range(8):
#      print("%6d,"%centermod(pinv*centermod(mont*pow(zeta,bitrev6(8*i+j),p),p),2**16),end='')
#    print()
#  print("    },")

  print("    .twist64 = {")
  for i in range(64//8):
    print("      ",end='')
    for j in range(8):
      print("%6d,"%centermod(mont*pow(zeta,8*i+j,p),p),end='')
    print()
  print("    },")
  print("    .twist32n = {")
  for i in range(32//8):
    print("      ",end='')
    for j in range(8):
      print("%6d,"%centermod(-mont*pow(zeta,2*(8*i+j),p),p),end='')
    print()
  print("    },")
  print("    .twist16 = {")
  for i in range(16//8):
    print("      ",end='')
    for j in range(8):
      print("%6d,"%centermod(mont*pow(zeta,4*(8*i+j),p),p),end='')
    print()
  print("    },")
  print("    .twist8n = {")
  print("      ",end='')
  for j in range(8):
    print("%6d,"%centermod(-mont*pow(zeta,8*j,p),p),end='')
  print()
  print("    },")
  print("    .twist4 = {")
  print("      ",end='')
  for j in range(4):
    print("%6d,"%centermod(mont*pow(zeta,16*j,p),p),end='')
  print()
  print("    },")
  print("    .twist2n = {")
  print("      ",end='')
  for j in range(2):
    print("%6d,"%centermod(-mont*pow(zeta,32*j,p),p),end='')
  print()
  print("    },")
  print("    .twist1 = {")
  print("      ",end='')
  for j in range(1):
    print("%6d,"%centermod(mont*pow(zeta,64*j,p),p),end='')
  print()
  print("    },")

  print("    .twist64_pinv = {")
  for i in range(64//8):
    print("      ",end='')
    for j in range(8):
      print("%6d,"%centermod(pinv*centermod(mont*pow(zeta,8*i+j,p),p),2**16),end='')
    print()
  print("    },")
  print("    .twist32n_pinv = {")
  for i in range(32//8):
    print("      ",end='')
    for j in range(8):
      print("%6d,"%centermod(pinv*centermod(-mont*pow(zeta,2*(8*i+j),p),p),2**16),end='')
    print()
  print("    },")
  print("    .twist16_pinv = {")
  for i in range(16//8):
    print("      ",end='')
    for j in range(8):
      print("%6d,"%centermod(pinv*centermod(mont*pow(zeta,4*(8*i+j),p),p),2**16),end='')
    print()
  print("    },")
  print("    .twist8n_pinv = {")
  print("      ",end='')
  for j in range(8):
    print("%6d,"%centermod(pinv*centermod(-mont*pow(zeta,8*j,p),p),2**16),end='')
  print()
  print("    },")
  print("    .twist4_pinv = {")
  print("      ",end='')
  for j in range(4):
    print("%6d,"%centermod(pinv*centermod(mont*pow(zeta,16*j,p),p),2**16),end='')
  print()
  print("    },")
  print("    .twist2n_pinv = {")
  print("      ",end='')
  for j in range(2):
    print("%6d,"%centermod(pinv*centermod(-mont*pow(zeta,32*j,p),p),2**16),end='')
  print()
  print("    },")
  print("    .twist1_pinv = {")
  print("      ",end='')
  for j in range(1):
    print("%6d,"%centermod(pinv*centermod(mont*pow(zeta,64*j,p),p),2**16),end='')
  print()
  print("    },")
  print("  },")

print("};")

#pmq = centermod(-P,q)
pmq = -P%q

print()
print("__attribute__((aligned(64)))")
print("const qdata modulus = {")
print("  .q = {{")
for i in range(math.ceil(nlimbs/8)):
  print("    ",end='')
  for j in range(8):
    if 8*i+j >= nlimbs: break
    print("%6d,"%((q >> (14*(8*i+j)))&0x3FFF),end='')
  print()
print("  }},")
print("  .pmq = {{")
for i in range(math.ceil(nlimbs/8)):
  print("    ",end='')
  for j in range(8):
    if 8*i+j < nlimbs-1:
      print("%6d,"%((pmq >> (14*(8*i+j)))&0x3FFF),end='')
    elif 8*i+j == nlimbs-1:
      print("%6d,"%(pmq >> (14*(8*i+j))),end='')
    else:
      break
  print()
print("  }},")
print("  .xvec = {")
for prime in primes[:nprimes]:
  p = prime[0]
  #x = centermod(P//p,q)
  x = P//p%q
  print("    {{")
  for i in range(math.ceil(nlimbs/8)):
    print("      ",end='')
    for j in range(8):
      if 8*i+j < nlimbs-1:
        print("%6d,"%((x >> (14*(8*i+j)))&0x3FFF),end='')
      elif 8*i+j == nlimbs-1:
        print("%6d,"%(x >> (14*(8*i+j))),end='')
      else:
        break
    print()
  print("    }},")
print("  },")
print("};")
