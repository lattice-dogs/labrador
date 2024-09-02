CC ?= /usr/bin/cc
CFLAGS += -std=c2x -Wall -Wextra -Wmissing-prototypes -Wredundant-decls \
  -Wshadow -Wpointer-arith -Wno-unused-function -flto=auto \
  -fwrapv -march=native -mtune=native -O3
RM = /bin/rm

SOURCES = pack.c greyhound.c dachshund.c chihuahua.c labrador.c \
  data.c jlproj.c polx.c poly.c polz.c sparsemat.c ntt.S invntt.S \
  aesctr.c fips202.c randombytes.c cpucycles.c
HEADERS = pack.h greyhound.h dachshund.h chihuahua.h labrador.h \
  data.h jlproj.h polx.h poly.h polz.h sparsemat.h fq.inc shuffle.inc \
  aesctr.h fips202.h randombytes.h malloc.h cpucycles.h

.PHONY: all

all: \
  test_aesctr \
  test_ntt \
  test_poly \
  test_polz \
  test_jlproj \
  test_chihuahua \
  test_dachshund \
  test_greyhound

%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

test_aesctr: test_aesctr.c aesctr.c aesctr.h randombytes.c randombytes.h cpucycles.c cpucycles.h
	$(CC) $(CFLAGS) test_aesctr.c aesctr.c randombytes.c cpucycles.c -o test_aesctr -lcrypto

test_ntt: test_ntt.c data.c data.h poly.c poly.h ntt.S invntt.S fq.inc shuffle.inc aesctr.c aesctr.h fips202.c fips202.h randombytes.c randombytes.h cpucycles.c cpucycles.h
	$(CC) $(CFLAGS) test_ntt.c data.c poly.c ntt.S invntt.S aesctr.c fips202.c randombytes.c cpucycles.c -o test_ntt -lm

test_poly: test_poly.c data.c data.h poly.c poly.h ntt.S invntt.S fq.inc shuffle.inc aesctr.c aesctr.h fips202.c fips202.h randombytes.c randombytes.h
	$(CC) $(CFLAGS) test_poly.c data.c poly.c ntt.S invntt.S aesctr.c fips202.c randombytes.c -o test_poly -lm

test_polz: test_polz.c data.c data.h polx.c polx.h poly.c poly.h polz.c polz.h ntt.S invntt.S fq.inc shuffle.inc aesctr.c aesctr.h fips202.c fips202.h randombytes.c randombytes.h cpucycles.c cpucycles.h
	$(CC) $(CFLAGS) test_polz.c data.c polx.c poly.c polz.c ntt.S invntt.S aesctr.c fips202.c randombytes.c cpucycles.c -o test_polz -lm -lgmp

test_jlproj: test_jlproj.c data.c data.h jlproj.c jlproj.h polx.c polx.h poly.c poly.h polz.c polz.h ntt.S invntt.S fq.inc shuffle.inc aesctr.c aesctr.h fips202.c fips202.h randombytes.c randombytes.h cpucycles.c cpucycles.h
	$(CC) $(CFLAGS) test_jlproj.c jlproj.c data.c polx.c poly.c polz.c ntt.S invntt.S aesctr.c fips202.c randombytes.c cpucycles.c -o test_jlproj -lm

test_chihuahua: test_chihuahua.c $(SOURCES) $(HEADERS)
	$(CC) $(CFLAGS) test_chihuahua.c $(SOURCES) -o $@ -lm

test_dachshund: test_dachshund.c $(SOURCES) $(HEADERS)
	$(CC) $(CFLAGS) test_dachshund.c $(SOURCES) -o $@ -lm

test_greyhound: test_greyhound.c $(SOURCES) $(HEADERS)
	$(CC) $(CFLAGS) test_greyhound.c $(SOURCES) -o $@ -lm

libdogs.so: $(SOURCES) $(HEADERS)
	$(CC) -shared -fPIC -fvisibility=hidden $(CFLAGS) -o $@ $(SOURCES)

clean:
	-$(RM) -rf *.o *.so
	-$(RM) -rf test_aesctr
	-$(RM) -rf test_ntt
	-$(RM) -rf test_poly
	-$(RM) -rf test_polz
	-$(RM) -rf test_jlproj
	-$(RM) -rf test_chihuahua
	-$(RM) -rf test_dachshund
	-$(RM) -rf test_greyhound
