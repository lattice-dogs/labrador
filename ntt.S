#include "data.h"
.include "shuffle.inc"
.include "fq.inc"

#define COMPRESS 0

.macro compress r0,r1
vpshufb		%ymm15,%ymm\r0,%ymm\r0
.if \r0 != \r1
vpshufb		%ymm15,%ymm\r1,%ymm\r1
vpunpcklqdq	%ymm\r1,%ymm\r0,%ymm\r0
.endif
vpermq		$0xD8,%ymm\r0,%ymm\r0
.endm

.macro update rln,rl0,rl1,rh0,rh1
vpaddw		%ymm\rh0,%ymm\rl0,%ymm\rln
vpsubw		%ymm\rh0,%ymm\rl0,%ymm\rh0
vpaddw		%ymm\rh1,%ymm\rl1,%ymm\rl0
vpsubw		%ymm\rh1,%ymm\rl1,%ymm\rh1
.endm

.macro levels0t5
/* twist */
#if COMPRESS
vmovdqa		 0*2(%rsi),%ymm1
vmovdqu		(_TWIST64_PINV+ 0)*2(%rdx),%ymm6
vmovdqu		(_TWIST64+ 0)*2(%rdx),%ymm7
fqmulprecomp	6,7,1,x=2
vmovdqa		16*2(%rsi),%ymm2
vmovdqu		(_TWIST64_PINV+16)*2(%rdx),%ymm8
vmovdqu		(_TWIST64+16)*2(%rdx),%ymm9
fqmulprecomp	8,9,2,x=3
vmovdqa		32*2(%rsi),%ymm3
vmovdqu		(_TWIST64_PINV+32)*2(%rdx),%ymm10
vmovdqu		(_TWIST64+32)*2(%rdx),%ymm11
fqmulprecomp	10,11,3,x=4
vmovdqa		48*2(%rsi),%ymm4
vmovdqu		(_TWIST64_PINV+48)*2(%rdx),%ymm12
vmovdqu		(_TWIST64+48)*2(%rdx),%ymm13
fqmulprecomp	12,13,4,x=5
#else
vmovdqa		 0*2(%rsi),%ymm1
fqmulprecomp2	(_TWIST64_PINV+ 0),(_TWIST64+ 0),1,x=2
vmovdqa		16*2(%rsi),%ymm2
fqmulprecomp2	(_TWIST64_PINV+16),(_TWIST64+16),2,x=3
vmovdqa		32*2(%rsi),%ymm3
fqmulprecomp2	(_TWIST64_PINV+32),(_TWIST64+32),3,x=4
vmovdqa		48*2(%rsi),%ymm4
fqmulprecomp2	(_TWIST64_PINV+48),(_TWIST64+48),4,x=5
#endif

/* level0 */
#if COMPRESS
//vmovdqa		_SHUFBIDX(%rip),%ymm15
movq		$0x0D0C090805040100,%rax
movq		%rax,%xmm15
vpbroadcastq	%xmm15,%ymm15
compress	6,8
compress	7,9
compress	10,12
compress	11,13
#endif

update		5,1,2,3,4

#if COMPRESS
fqmulprecomp	6,7,3,x=2
fqmulprecomp	10,11,4,x=2
#else
fqmulprecomp2	(_TWIST32N_PINV+ 0),(_TWIST32N+ 0),3,x=2,neg=1
fqmulprecomp2	(_TWIST32N_PINV+16),(_TWIST32N+16),4,x=2,neg=1
#endif

//5,1,3,4

/* level1 */
#if COMPRESS
compress	6,10
compress	7,11
#else
vmovdqu		_TWIST16_PINV*2(%rdx),%ymm6
vmovdqu		_TWIST16*2(%rdx),%ymm7
#endif
vpbroadcastw	%xmm6,%ymm8
vpbroadcastw	%xmm7,%ymm9

update		2,5,3,1,4
fqmulprecomp	6,7,1,x=3
fqmulprecomp	8,9,2,x=3
fqmulprecomp	6,7,4,x=3
fqmulprecomp	8,9,5,x=3

//2,1,5,4

shuffle8	2,5,3,5
shuffle8	1,4,2,4

//3,5,2,4

/* level2 */
#if COMPRESS
compress	6,6
compress	7,7
#else
vbroadcasti128	_TWIST8N_PINV*2(%rdx),%ymm6
vbroadcasti128	_TWIST8N*2(%rdx),%ymm7
#endif

update		1,3,2,5,4
#if COMPRESS
fqmulprecomp	6,7,5,x=2
fqmulprecomp	6,7,4,x=2
#else
fqmulprecomp	6,7,5,x=2,neg=1
fqmulprecomp	6,7,4,x=2,neg=1
#endif

//1,5,3,4

shuffle4	1,3,2,3
shuffle4	5,4,1,4

//2,3,1,4

/* level3 */
#if COMPRESS
vpshufb		%ymm15,%ymm6,%ymm6
vpshufb		%ymm15,%ymm7,%ymm7
#else
vpbroadcastq	_TWIST4_PINV*2(%rdx),%ymm6
vpbroadcastq	_TWIST4*2(%rdx),%ymm7
#endif

update		5,2,1,3,4
fqmulprecomp	6,7,3,x=1
fqmulprecomp	8,9,5,x=1
fqmulprecomp	6,7,4,x=1
fqmulprecomp	8,9,2,x=1

//5,3,2,4

shuffle2	5,2,1,2
shuffle2	3,4,5,4

//1,2,5,4

/* level4 */
#if COMPRESS
vpshufb		%ymm15,%ymm6,%ymm6
vpshufb		%ymm15,%ymm7,%ymm7
#else
vpbroadcastd	_TWIST2N_PINV*2(%rdx),%ymm6
vpbroadcastd	_TWIST2N*2(%rdx),%ymm7
#endif

update		3,1,5,2,4
#if COMPRESS
fqmulprecomp	6,7,2,x=5
fqmulprecomp	6,7,4,x=5
#else
fqmulprecomp	6,7,2,x=5,neg=1
fqmulprecomp	6,7,4,x=5,neg=1
#endif

//3,2,1,4

shuffle1	3,1,5,1
shuffle1	2,4,3,4

//5,1,3,4

/* level5 */
vpbroadcastw	_PINV*2(%rdx),%ymm6
vpbroadcastw	_F*2(%rdx),%ymm7
vpmullw		%ymm6,%ymm7,%ymm6

update		2,5,3,1,4
fqmulprecomp	6,7,2,x=3
fqmulprecomp	6,7,1,x=3
fqmulprecomp	6,7,5,x=3
fqmulprecomp	6,7,4,x=3

//2,1,5,4

vmovdqa		%ymm2, 0*2(%rdi)
vmovdqa		%ymm1,16*2(%rdi)
vmovdqa		%ymm5,32*2(%rdi)
vmovdqa		%ymm4,48*2(%rdi)
.endm

.text
.global poly_ntt
poly_ntt:
vpbroadcastw	_P*2(%rdx),%ymm0

levels0t5

ret

.section .note.GNU-stack,"",@progbits
