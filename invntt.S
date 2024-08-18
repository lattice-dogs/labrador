#include "data.h"
.include "shuffle.inc"
.include "fq.inc"

.macro update rln,rl0,rl1,rh0,rh1
vpaddw		%ymm\rh0,%ymm\rl0,%ymm\rln
vpsubw		%ymm\rh0,%ymm\rl0,%ymm\rh0
vpaddw		%ymm\rh1,%ymm\rl1,%ymm\rl0
vpsubw		%ymm\rh1,%ymm\rl1,%ymm\rh1
.endm

.macro ilevels0t5
/* level0 */
vmovdqa		 0*2(%rsi),%ymm1
vmovdqa		16*2(%rsi),%ymm2
vmovdqa		32*2(%rsi),%ymm3
vmovdqa		48*2(%rsi),%ymm4

update		5,1,3,2,4

//[5,2],[1,4]

/* level1 */
vpbroadcastd	(_TWIST2N_PINV+1)*2(%rdx),%ymm6
vpbroadcastd	(_TWIST2N+1)*2(%rdx),%ymm7
vpshufb		%ymm15,%ymm6,%ymm6
vpshufb		%ymm15,%ymm7,%ymm7

shuffle1	5,2,3,2
shuffle1	1,4,5,4

//3,5,2,4

fqmulprecomp	6,7,5,x=1,neg=0
fqmulprecomp	6,7,4,x=1,neg=0

update		1,3,2,5,4

//[1,5],[3,4]

/* level2 */
vpbroadcastq	(_TWIST4_PINV+1)*2(%rdx),%ymm6
vpbroadcastq	(_TWIST4+1)*2(%rdx),%ymm7
vpshufb		%ymm15,%ymm6,%ymm6
vpshufb		%ymm15,%ymm7,%ymm7

shuffle2	1,5,2,5
shuffle2	3,4,1,4

//2,1,5,4

fqmulprecomp	6,7,1,x=3,neg=1
fqmulprecomp	6,7,4,x=3,neg=1

update		3,2,5,1,4

//[3,1],[2,4]

/* level3 */
vbroadcasti128	(_TWIST8N_PINV+1)*2(%rdx),%ymm6
vbroadcasti128	(_TWIST8N+1)*2(%rdx),%ymm7
vpshufb		%ymm15,%ymm6,%ymm6
vpshufb		%ymm15,%ymm7,%ymm7
vpbroadcastw	(_TWIST16_PINV)*2(%rdx),%ymm8
vpbroadcastw	(_TWIST16)*2(%rdx),%ymm9

shuffle4	3,1,5,1
shuffle4	2,4,3,4

//5,3,1,4

fqmulprecomp	6,7,3,x=2,neg=0
fqmulprecomp	6,7,4,x=2,neg=0
fqmulprecomp	8,9,5,x=2
fqmulprecomp	8,9,1,x=2

update		2,5,1,3,4

//[2,3],[5,4]

/* level4 */
vmovdqu		(_TWIST16_PINV+1)*2(%rdx),%ymm6
vmovdqu		(_TWIST16+1)*2(%rdx),%ymm7
vpshufb		%ymm15,%ymm6,%ymm6
vpshufb		%ymm15,%ymm7,%ymm7
vperm2i128	$0x01,%ymm6,%ymm6,%ymm6
vperm2i128	$0x01,%ymm7,%ymm7,%ymm7

shuffle8	2,3,1,3
shuffle8	5,4,2,4

//1,2,3,4

fqmulprecomp	6,7,2,x=5,neg=1
fqmulprecomp	6,7,4,x=5,neg=1

update		5,1,3,2,4

//[5,2],[1,4]

/* level5 */
vmovdqu		(_TWIST32N_PINV+17)*2(%rdx),%ymm6
vmovdqu		(_TWIST32N+17)*2(%rdx),%ymm7
vmovdqu		(_TWIST32N_PINV+ 1)*2(%rdx),%ymm8
vmovdqu		(_TWIST32N+ 1)*2(%rdx),%ymm9
vpshufb		%ymm15,%ymm6,%ymm6
vpshufb		%ymm15,%ymm7,%ymm7
vpshufb		%ymm15,%ymm8,%ymm8
vpshufb		%ymm15,%ymm9,%ymm9
vperm2i128	$0x01,%ymm6,%ymm6,%ymm6
vperm2i128	$0x01,%ymm7,%ymm7,%ymm7
vperm2i128	$0x01,%ymm8,%ymm8,%ymm8
vperm2i128	$0x01,%ymm9,%ymm9,%ymm9

fqmulprecomp	6,7,1,x=3,neg=0
fqmulprecomp	8,9,4,x=3,neg=0

update		3,5,2,1,4

//[3,5,1,4]

/* Twist */
vmovdqu		(_TWIST64_PINV+49)*2(%rdx),%ymm6
vmovdqu		(_TWIST64+49)*2(%rdx),%ymm7
vmovdqu		(_TWIST64_PINV+33)*2(%rdx),%ymm8
vmovdqu		(_TWIST64+33)*2(%rdx),%ymm9
vmovdqu		(_TWIST64_PINV+17)*2(%rdx),%ymm10
vmovdqu		(_TWIST64+17)*2(%rdx),%ymm11
vmovdqu		(_TWIST64_PINV+ 1)*2(%rdx),%ymm12
vmovdqu		(_TWIST64+ 1)*2(%rdx),%ymm13
vpshufb		%ymm15,%ymm6,%ymm6
vpshufb		%ymm15,%ymm7,%ymm7
vpshufb		%ymm15,%ymm8,%ymm8
vpshufb		%ymm15,%ymm9,%ymm9
vpshufb		%ymm15,%ymm10,%ymm10
vpshufb		%ymm15,%ymm11,%ymm11
vpshufb		%ymm15,%ymm12,%ymm12
vpshufb		%ymm15,%ymm13,%ymm13
vperm2i128	$0x01,%ymm6,%ymm6,%ymm6
vperm2i128	$0x01,%ymm7,%ymm7,%ymm7
vperm2i128	$0x01,%ymm8,%ymm8,%ymm8
vperm2i128	$0x01,%ymm9,%ymm9,%ymm9
vperm2i128	$0x01,%ymm10,%ymm10,%ymm10
vperm2i128	$0x01,%ymm11,%ymm11,%ymm11
vperm2i128	$0x01,%ymm12,%ymm12,%ymm12
vperm2i128	$0x01,%ymm13,%ymm13,%ymm13

fqmulprecomp	6,7,3,x=2,neg=1
fqmulprecomp	8,9,5,x=2,neg=1
fqmulprecomp	10,11,1,x=2,neg=1
fqmulprecomp	12,13,4,x=2,neg=1

vmovdqa		%ymm3, 0*2(%rdi)
vmovdqa		%ymm5,16*2(%rdi)
vmovdqa		%ymm1,32*2(%rdi)
vmovdqa		%ymm4,48*2(%rdi)
.endm

.text
.global poly_invntt
poly_invntt:
vpbroadcastw	_P*2(%rdx),%ymm0

movq		$0x09080B0A0D0C0F0E,%rax
movq		%rax,%xmm15
movq		$0x0100030205040706,%rax
pinsrq		$0x1,%rax,%xmm15
vinserti128	$0x1,%xmm15,%ymm15,%ymm15

ilevels0t5

ret

.section .note.GNU-stack,"",@progbits
