// Used as an X-macro to build various Fitch algorithm variants.
// This file should not be #included directly, rather you should #include "generate_all_fitch_b2_w4_fastc.inc",
// which will build appropriately-named declarations (if PARAM_HEADER is #defined) or definitions
// (if not) for all combinations of STOPEARLY and FASTWEIGHT1FITCH.
// Also the calling code must #include "switches.h" so that FASTCALL is defined appropriately.

#define PASTE2(x, y) x ## y
#define PASTE(x, y) PASTE2(x, y)

// Make sure we're not stamping on anyone's toes
#if defined(TEMP_SUFFIX1) || defined(TEMP_SUFFIX2) || defined(TEMP_SUFFIX3)
#error "At least one of the macros TEMP_SUFFIX1, TEMP_SUFFIX2 and TEMP_SUFFIX3 is already in use!"
#endif	// defined(TEMP_SUFFIX1) || defined(TEMP_SUFFIX2) || defined(TEMP_SUFFIX3)

#define TEMP_SUFFIX1 _b2_w16_sse2asm

#ifdef PARAM_STOPEARLY
#define TEMP_SUFFIX2 PASTE(_stopearly, TEMP_SUFFIX1)
//#pragma message("PARAM_STOPEARLY")
#else	// not PARAM_STOPEARLY
#define TEMP_SUFFIX2 TEMP_SUFFIX1
//#pragma message("no PARAM_STOPEARLY")
#endif	// not PARAM_STOPEARLY

#ifdef PARAM_FASTWEIGHT1FITCH
#define TEMP_SUFFIX3 PASTE(_fw1, TEMP_SUFFIX2)
//#pragma message("PARAM_FASTWEIGHT1FITCH")
#else	// not PARAM_FASTWEIGHT1FITCH
#define TEMP_SUFFIX3 TEMP_SUFFIX2
//#pragma message("no PARAM_FASTWEIGHT1FITCH")
#endif	// not PARAM_FASTWEIGHT1FITCH


///////////////////////////////////////////////////////////////////////////////
// SSE2ASM routines here.
///////////////////////////////////////////////////////////////////////////////


// WTJW 16/2/2005: I think all combinations of STOPEARLY, FASTWEIGHT1FITCH and
// WEIGHT1MULTIPLY work, both on Windows under MSVC++ and on Linux under gcc.
// Regarding use of esi and edi regs: MSVC++ promises not to allocate these to
// local vars when they are used by asm code (which is fine since there are no
// local vars here) and gcc is notified that we use them by the (script-created)
// clobber list, and works around them as need be.  Thus there is no need to
// PUSH and POP these regs.
unsigned FASTCALL PASTE(FitchScore, TEMP_SUFFIX3)(void *s1, void *s2, unsigned *weights, unsigned nBlocks
#ifdef PARAM_STOPEARLY
	, unsigned maxScore
#endif	// PARAM_STOPEARLY
)
#ifdef PARAM_HEADER
;
#else	// note PARAM_HEADER
{
	SSE2ALIGN static unsigned xmm7val[4] = { 0x00000000, 0x0000000F, 0x000000FF, 0x00000FFF };		// Mask for ">" comparison
	SSE2ALIGN static unsigned xmm6val[4] = { 0x0000000F, 0x000000F0, 0x00000F00, 0x0000F000 };		// Mask for ANDing
#ifdef PARAM_FASTWEIGHT1FITCH
	SSE2ALIGN static unsigned xmm7valW1[4] = { 0x77777777, 0x77777777, 0x77777777, 0x77777777 };
//	SSE2ALIGN static unsigned xmm6valW1[4] = { 0x11111111, 0x11111111, 0x11111111, 0x11111111 };
	SSE2ALIGN static unsigned xmm6valW1[4] = { 0x88888888, 0x88888888, 0x88888888, 0x88888888 };
//	SSE2ALIGN static unsigned xmm3valW1[4] = { 0x0000000F, 0x0000000F, 0x0000000F, 0x0000000F };
//	SSE2ALIGN static unsigned xmm3valW1[4] = { 0x0000000F, 0x00000000, 0x0000000F, 0x00000000 };
//	SSE2ALIGN static unsigned xmm3valW1[4] = { 0x00000078, 0x00000078, 0x00000078, 0x00000078 };
	SSE2ALIGN static unsigned xmm3valW1[4] = { 0x00000078, 0x00000000, 0x00000078, 0x00000000 };
#endif	// PARAM_FASTWEIGHT1FITCH
//	SSE2ALIGN static unsigned temp[4] = { 0x00000000, 0x00000000, 0x00000000, 0x00000000 };		//DEBUG
//	static unsigned temp = 0;		//DEBUG
#ifdef PARAM_STOPEARLY
#ifdef PARAM_FASTWEIGHT1FITCH
	SSE2ALIGN static unsigned xmm5valW1[4] = { 0xFFFFFFFF, 0x00000000, 0x00000000, 0x00000000 };
#endif	// PARAM_FASTWEIGHT1FITCH
//	SSE2ALIGN static unsigned xmmEvenDwordMask[4] = { 0xFFFFFFFF, 0x00000000, 0xFFFFFFFF, 0x00000000 };
#endif	// PARAM_STOPEARLY
	unsigned thisScore = 0;
	
	INCMEASURE(scoreFitchCount);
	
//	temp = 0;		//DEBUG
	
	__asm {
#ifdef PARAM_STOPEARLY
//		push	esi
//		push	edi
		xor		esi,esi			// By the time end1 is reached, this must contain the total.
		mov		edi,[maxScore]	// {{GCC I "D" (maxScore)}} {{GCC DEL}}
//		movdqa	xmm5,[xmm5val]	// {{GCC I [xmm5val] "m" (*xmm5val)}} {{GCC INS "movdqa %[xmm5val],%%xmm5\n\t"}} {{GCC DEL}}
#endif	// not PARAM_STOPEARLY
		
		movdqa	xmm7,[xmm7val]	// {{GCC I [xmm7val] "m" (*xmm7val)}} {{GCC INS "movdqa %[xmm7val],%%xmm7\n\t"}} {{GCC DEL}}
		movdqa	xmm6,[xmm6val]	// {{GCC I [xmm6val] "m" (*xmm6val)}} {{GCC INS "movdqa %[xmm6val],%%xmm6\n\t"}} {{GCC DEL}}
		pxor	xmm4,xmm4		// Holds total
//		pxor	xmm5,xmm5		// Another total
		
		mov		ecx,[nBlocks]	// {{GCC I "c" (nBlocks)}} {{GCC DEL}}
		add		ecx,ecx
		add		ecx,ecx
		add		ecx,ecx
		add		ecx,ecx			// ecx now counts bytes
		mov		eax,[s1]		// {{GCC I "a" (s1)}} {{GCC DEL}}
		mov		ebx,[s2]		// {{GCC I "b" (s2)}} {{GCC DEL}}
		mov		edx,[weights]	// {{GCC I "d" (weights)}} {{GCC DEL}}
		add		eax,ecx
		add		ebx,ecx			// NB: We increment edx separately from ecx.
		neg		ecx				// Now [esi+ecx] etc. point to the 1st element
		jz		end1			// Exit at this point if nothing to do
		
	top1:
#ifdef PARAM_FASTWEIGHT1FITCH
		cmp		[edx],1			// Stupid intel2gas converts this instruction wrongly, so do it manually: {{GCC DEL}} {{GCC INS "cmpl $1,(%%edx)\n\t"}}
		jne		notweight1
		
		// We are dealing with weight-1 columns from hereonin
		
		// Setup
#ifdef PARAM_STOPEARLY
//		push	esi
//		push	edi
		// Since the PSADBW instruction creates only two qword values (rather
		// than 4 dword values), let's combine the score so far into 2 qwords.
		pxor	xmm4,xmm4
		movd	xmm4,esi
		movdqa	xmm5,[xmm5valW1]	// {{GCC I [xmm5valW1] "m" (*xmm5valW1)}} {{GCC INS "movdqa %[xmm5valW1],%%xmm5\n\t"}} {{GCC DEL}}
		shl		edi,3			// Compensate for left shift that score gets
#endif	// PARAM_STOPEARLY
		
		pslld	xmm4,3			// We will shift right by this amount at the end, so compensate now
		movdqa	xmm7,[xmm7valW1]	// {{GCC I [xmm7valW1] "m" (*xmm7valW1)}} {{GCC INS "movdqa %[xmm7valW1],%%xmm7\n\t"}} {{GCC DEL}}
		movdqa	xmm6,[xmm6valW1]	// {{GCC I [xmm6valW1] "m" (*xmm6valW1)}} {{GCC INS "movdqa %[xmm6valW1],%%xmm6\n\t"}} {{GCC DEL}}
		movdqa	xmm3,[xmm3valW1]	// {{GCC I [xmm3valW1] "m" (*xmm3valW1)}} {{GCC INS "movdqa %[xmm3valW1],%%xmm3\n\t"}} {{GCC DEL}}
		
	weight1:
		movdqa	xmm0,[eax+ecx]
//		movdqa	xmm1,[ebx+ecx]
//		pand	xmm0,xmm1		// Put intersection in xmm0
		pand	xmm0,[ebx+ecx]		// Put intersection in xmm0
		movdqa	xmm2,xmm0
		
		pand	xmm0,xmm7
		paddb	xmm0,xmm7
		por		xmm0,xmm2
//		psrld	xmm0,3
		pand	xmm0,xmm6		// Now the low bit is OFF in each nibble iff there were no bits on in that nibble
#ifndef PARAM_STOPEARLY
		add		ecx,16			// Moved ahead to use pipeline gaps (MMX ops don't affect flags)
#endif	// not PARAM_STOPEARLY
		
		// For each of low QWORD and high QWORD: Sum (1-x) for each of the 8 low
		// nibbles into low nibble of QWORD and then do the same for the 8 high
		// nibbles into high nibble of low byte of QWORD!
		psadbw	xmm0,xmm6
		movdqa	xmm1,xmm0
		pand	xmm1,xmm3		// Now we just have a weight in the low byte of each QWORD in xmm0.
		psrld	xmm0,4
		pand	xmm0,xmm3
		
#ifdef PARAM_STOPEARLY
		paddd	xmm1,xmm0
		paddd	xmm4,xmm1
		pshufd	xmm0,xmm1,10101010b	// Place dword 2 in dword 0 of xmm0
		paddd	xmm4,xmm0
		
		movd	esi,xmm4
		cmp		esi,edi
		ja		shiftback
		
		pand	xmm4,xmm5		// Set high 3 dwords to zero so total is in bottom dword
		add		ecx,16
#else	// not PARAM_STOPEARLY
		paddd	xmm4,xmm1
		paddd	xmm4,xmm0
#endif	// not PARAM_STOPEARLY

//#ifdef WEIGHT1MULTIPLY
//		pmaddwd	xmm0,xmm6		// Voila!  Sum of each WORD's bits is in nibble 3 of each DWORD.
//		psrld	xmm0,12
//		pand	xmm0,xmm3
//#else	// not WEIGHT1MULTIPLY
//		movdqa	xmm1,xmm0
//		psrld	xmm0,4
//		paddb	xmm0,xmm1		// Adjacent nibbles have now been summed
//		movdqa	xmm1,xmm0
//		psrld	xmm0,8
//		paddb	xmm0,xmm1		// Adjacent bytes have now been summed
//		movdqa	xmm1,xmm0
//		psrld	xmm0,16
//		paddb	xmm0,xmm1		// Adjacent words have now been summed
//		pand	xmm0,xmm3
//#endif	// not WEIGHT1MULTIPLY
//		paddd	xmm4,xmm0			// Add to total #1
		
		jnz		weight1
		
	shiftback:
#ifdef PARAM_STOPEARLY
		shr		esi,3
#else	// not PARAM_STOPEARLY
		psrld	xmm4,3
#endif	// not PARAM_STOPEARLY
		jmp		end1
		
	notweight1:
#endif	// PARAM_FASTWEIGHT1FITCH
		// Strategy: We get two independent "threads" of data-dependency:
		// The first 8 nibbles, and the 2nd 8 nibbles.  Hopefully, P4
		// pipelining will enable these "threads" to be efficiently "mingled."
		movdqa	xmm0,[eax+ecx]
		movdqa	xmm1,[ebx+ecx]
		pand	xmm0,xmm1		// Put intersection in xmm0
		
		pshufd	xmm1,xmm0,0		// Copy DWORD 0 to all four DWORDS in xmm1
		pshufd	xmm2,xmm0,055h	// Copy DWORD 1 to all four DWORDS in xmm2
		pshufd	xmm3,xmm0,0AAh	// Copy DWORD 2 to all four DWORDS in xmm3
		pshufd	xmm0,xmm0,0FFh	// Copy DWORD 3 to all four DWORDS in xmm0 itself
		
		movdqa	xmm5,xmm1
		psrld	xmm5,16			// In preparation for nibbles 4, 5, 6 and 7
		
		// Process nibbles 0, 1 2 and 3
		pand	xmm1,xmm6
		pcmpgtd	xmm1,xmm7		// DWORD 0 tests low nibble, DWORD 1 tests next nibble, ..., DWORD 3 tests high nibble in 16-bit word
		pandn	xmm1,[edx]
		paddd	xmm4,xmm1		// Accumulate total
		
		// Process nibbles 4, 5, 6 and 7
		pand	xmm5,xmm6
		pcmpgtd	xmm5,xmm7		// DWORD 0 tests low nibble, DWORD 1 tests next nibble, ..., DWORD 3 tests high nibble in 16-bit word
		pandn	xmm5,[edx+16]
		paddd	xmm4,xmm5		// Accumulate total
		
		movdqa	xmm5,xmm2
		psrld	xmm5,16			// In preparation for nibbles 12, 13, 14 and 15
		
		// Process nibbles 8, 9, 10 and 11
		pand	xmm2,xmm6
		pcmpgtd	xmm2,xmm7		// DWORD 0 tests low nibble, DWORD 1 tests next nibble, ..., DWORD 3 tests high nibble in 16-bit word
		pandn	xmm2,[edx+32]
		paddd	xmm4,xmm2		// Accumulate total
		
		// Process nibbles 12, 13, 14 and 15
		pand	xmm5,xmm6
		pcmpgtd	xmm5,xmm7		// DWORD 0 tests low nibble, DWORD 1 tests next nibble, ..., DWORD 3 tests high nibble in 16-bit word
		pandn	xmm5,[edx+48]
		paddd	xmm4,xmm5		// Accumulate total
		
		movdqa	xmm5,xmm3
		psrld	xmm5,16			// In preparation for nibbles 20, 21, 22 and 23
		
		// Process nibbles 16, 17, 18 and 19
		pand	xmm3,xmm6
		pcmpgtd	xmm3,xmm7		// DWORD 0 tests low nibble, DWORD 1 tests next nibble, ..., DWORD 3 tests high nibble in 16-bit word
		pandn	xmm3,[edx+64]
		paddd	xmm4,xmm3		// Accumulate total
		
		// Process nibbles 20, 21, 22 and 23
		pand	xmm5,xmm6
		pcmpgtd	xmm5,xmm7		// DWORD 0 tests low nibble, DWORD 1 tests next nibble, ..., DWORD 3 tests high nibble in 16-bit word
		pandn	xmm5,[edx+80]
		paddd	xmm4,xmm5		// Accumulate total
		
		movdqa	xmm5,xmm0
		psrld	xmm5,16			// In preparation for nibbles 4, 5, 6 and 7
		
		// Process nibbles 24, 25, 26 and 27
		pand	xmm0,xmm6
		pcmpgtd	xmm0,xmm7		// DWORD 0 tests low nibble, DWORD 1 tests next nibble, ..., DWORD 3 tests high nibble in 16-bit word
		pandn	xmm0,[edx+96]
		paddd	xmm4,xmm0		// Accumulate total
		
		// Process nibbles 28, 29, 30 and 31
		pand	xmm5,xmm6
		pcmpgtd	xmm5,xmm7		// DWORD 0 tests low nibble, DWORD 1 tests next nibble, ..., DWORD 3 tests high nibble in 16-bit word
		pandn	xmm5,[edx+112]
		paddd	xmm4,xmm5		// Accumulate total
		
#ifdef PARAM_STOPEARLY
		movdqa	xmm5,xmm4		// Use xmm5 as a scratch reg for computing total
		pshufd	xmm0,xmm4,031h
		paddd	xmm5,xmm0		// Add adjacent pairs of DWORDs
		pshufd	xmm0,xmm5,2		// Add those sums to get a final total in the bottom DWORD
		paddd	xmm5,xmm0
		movd	esi,xmm5
		cmp		esi,edi
		ja		end1
#endif	// PARAM_STOPEARLY
		
		sub		edx,-128		// Just like "add edx,128" but without the fat :) (-128 can fit in 1 sign-extended byte)
		add		ecx,16
		jnz		top1
		
	end1:
//		paddd	xmm4,xmm5
		
//		pshufd	xmm0,xmm4,1
//		pshufd	xmm1,xmm4,2
//		pshufd	xmm2,xmm4,3
//		paddd	xmm4,xmm0
//		paddd	xmm4,xmm1
//		paddd	xmm4,xmm2
		
#ifdef PARAM_STOPEARLY
		mov		[thisScore],esi	// {{GCC O "=m" (thisScore)}} {{GCC DEL}} {{GCC INS "mov %%esi,%0\n\t"}}
#else	// not PARAM_STOPEARLY
		pshufd	xmm0,xmm4,031h
		paddd	xmm4,xmm0		// Add adjacent pairs of DWORDs
		pshufd	xmm0,xmm4,2		// Add those sums to get a final total in the bottom DWORD
		paddd	xmm4,xmm0
		movd	[thisScore],xmm4	// {{GCC O "=m" (thisScore)}} {{GCC DEL}} {{GCC INS "movd %%xmm4,%0\n\t"}}
#endif	// not PARAM_STOPEARLY
		
//		movq	mm0,mm4
//		punpckhdq	mm4,mm4		// Copy high DWORD to low DWORD
//		paddd	mm4,mm0			// Final total now in low DWORD of mm4
//		movd	[thisScore],mm4	// {{GCC O "=m" (thisScore)}} {{GCC DEL}} {{GCC INS "movd %%mm4,%0\n\t"}}
	}
	
	return thisScore;
}
#endif	// not PARAM_HEADER

// Clean up
#undef TEMP_SUFFIX1
#undef TEMP_SUFFIX2
#undef TEMP_SUFFIX3
