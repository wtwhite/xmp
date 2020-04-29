#include <assert.h>
#include "common.h"
#include <dbgprint.h>
#include "seq_b2_sse2asm.h"
#include "measure.h"
//#include "yanbader.h"

//#define SSE2ALIGN
#ifdef _MSC_VER
#define SSE2ALIGN __declspec(align(16))
#else	// not _MSC_VER
#define SSE2ALIGN __attribute__((aligned(16)))
//#define SSE2ALIGN "FOONLY"
#endif	// not _MSC_VER

///////////////////////////////////////////////////////////////////////////////
// SSE2ASM routines here.
///////////////////////////////////////////////////////////////////////////////

SSE2ALIGN static unsigned lowNibblesOn[4] = { 0x0F0F0F0F, 0x0F0F0F0F, 0x0F0F0F0F, 0x0F0F0F0F };

void FASTCALL FitchBases_b2_w16_sse2asm(void *s1, void *s2, void *dest, unsigned nBlocks) {
	INCMEASURE(fitchBasesCount);
	
//	assert((len & 3) == 0);
//	assert((((unsigned) xmm7val) & 15) == 0);
//	DBGPRINT1("START\n");
	
	__asm {
		// Using __fastcall, so we expect s1 in ecx and s2 in edx, rest on the stack.
		mov		eax,[nBlocks]	// {{GCC I "a" (nBlocks)}} {{GCC DEL}}
		add		eax,eax
		add		eax,eax
		add		eax,eax
		add		eax,eax			// eax now counts bytes
		jz		end1			// Exit at this point if nothing to do
		
		neg		eax
		mov		ebx,[dest]		// {{GCC I "b" (dest)}} {{GCC DEL}}
		sub		ecx,eax			// Now [ecx+eax] points to the 1st element of s1
		sub		edx,eax			// Now [edx+eax] points to the 1st element of s2
		sub		ebx,eax			// Now [ebx+eax] points to the 1st element of dest

		movdqa	xmm7,[lowNibblesOn]	// {{GCC I [lowNibblesOn] "m" (*lowNibblesOn)}} {{GCC INS "movdqa %[lowNibblesOn],%%xmm7\n\t"}} {{GCC DEL}}
		pxor	xmm5,xmm5		// A constant 0
		
	top1:
		movdqa	xmm0,[eax+ecx]
		movdqa	xmm1,[eax+edx]
		movdqa	xmm2,xmm0
		add		eax,16			// Moved ahead to fill pipeline gaps; OK since SSE2 operations don't affect flags

		pand	xmm0,xmm1		// Put intersection in xmm0
		por		xmm1,xmm2		// Put union in xmm1

		movdqa	xmm3,xmm0		// Save intersection out to xmm3
		movdqa	xmm2,xmm0		// Save intersection out to xmm2

		pand	xmm0,xmm7
		pcmpeqb	xmm3,xmm0		// Each byte is all 1s iff hi nibble had NO bits on (see why? :))
		pcmpeqb	xmm0,xmm5		// Each byte is all 1s iff low nibble had NO bits on
		por		xmm3,xmm7		// Turn low nibble bits on for combining later
		pandn	xmm0,xmm7		// Low-nibble bits off iff low nibble had NO bits on; hi-nibble bits OFF

		pandn	xmm0,xmm3		// Final union mask: Each nibble on iff no bit in that nibble was on
		pand	xmm0,xmm1		// Needed unions in xmm0
		por		xmm0,xmm2		// Final result in xmm0

		movdqa	[eax+ebx-16],xmm0
		jnz		top1
		
	end1:
	}
//	DBGPRINT1("END\n");
}

// Same as FitchBases_b2_w16_sse2asm(), except it checks whether the just-created sequence is the
// same as a pre-existing sequence at the same time.
unsigned FASTCALL FitchBasesAndCompare_b2_w16_sse2asm(void *s1, void *s2, void *dest, void *prev, unsigned nBlocks) {
	INCMEASURE(fitchBasesCount);
	
//	assert((len & 3) == 0);
//	assert((((unsigned) xmm7val) & 15) == 0);
//	DBGPRINT1("START\n");
	
	__asm {
		// Using __fastcall, so we expect s1 in ecx and s2 in edx, rest on the stack.
		mov		eax,[nBlocks]	// {{GCC I "a" (nBlocks)}} {{GCC DEL}}
		add		eax,eax
		add		eax,eax
		add		eax,eax
		add		eax,eax			// eax now counts bytes

		pxor	xmm5,xmm5		// A constant 0
		pxor	xmm4,xmm4		// Will hold 0 at the end iff there were no differences
		jz		end1			// Exit at this point if nothing to do

		neg		eax
		mov		esi,[prev]		// {{GCC I "S" (prev)}} {{GCC DEL}}
		mov		ebx,[dest]		// {{GCC I "b" (dest)}} {{GCC DEL}}
		sub		ecx,eax			// Now [ecx+eax] points to the 1st element of s1
		sub		edx,eax			// Now [edx+eax] points to the 1st element of s2
		sub		esi,eax			// Now [esi+eax] points to the 1st element of prev
		sub		ebx,eax			// Now [ebx+eax] points to the 1st element of dest
		
		movdqa	xmm7,[lowNibblesOn]	// {{GCC I [lowNibblesOn] "m" (*lowNibblesOn)}} {{GCC INS "movdqa %[lowNibblesOn],%%xmm7\n\t"}} {{GCC DEL}}

	top1:
		movdqa	xmm0,[eax+ecx]
		movdqa	xmm1,[eax+edx]
		movdqa	xmm6,[eax+esi]	// Moved ahead (for no particular reason...)
		add		eax,16			// Moved ahead to fill pipeline gaps; OK since SSE2 operations don't affect flags
		
		movdqa	xmm2,xmm0
		pand	xmm0,xmm1		// Put intersection in xmm0
		por		xmm1,xmm2		// Put union in xmm1
		
		movdqa	xmm3,xmm0		// Save intersection out to xmm3
		movdqa	xmm2,xmm0		// Save intersection out to xmm2

		pand	xmm0,xmm7
		pcmpeqb	xmm3,xmm0		// Each byte is all 1s iff hi nibble had NO bits on (see why? :))
		pcmpeqb	xmm0,xmm5		// Each byte is all 1s iff low nibble had NO bits on
		por		xmm3,xmm7		// Turn low nibble bits on for combining later
		pandn	xmm0,xmm7		// Low-nibble bits off iff low nibble had NO bits on; hi-nibble bits OFF

		pandn	xmm0,xmm3		// Final union mask: Each nibble on iff no bit in that nibble was on
		pand	xmm0,xmm1		// Needed unions in xmm0
		por		xmm0,xmm2		// Final result in xmm0

		movdqa	[eax+ebx-16],xmm0
		pxor	xmm6,xmm0		// xmm6 will be 0 iff all bytes are identical
		por		xmm4,xmm6
		jnz		top1
		
	end1:
		pcmpeqb	xmm4,xmm5		// Turns on all bits in bytes wherever the byte was 0
		pmovmskb	eax,xmm4	// High latency on P4 sadly, but faster than movd!!
		sub		eax,65535		// eax == 0 iff xmm4 was 0 (i.e. if there were no differences)
	}

	// No explicit return, since return value is in eax.
}

unsigned CompareSeqs_b2_w16_sse2asm(void *s1, void *s2, unsigned nBlocks) {
	unsigned retVal;

	__asm {
		mov		ecx,[nBlocks]	// {{GCC I "c" (nBlocks)}} {{GCC DEL}}
		xor		edx,edx
		add		ecx,ecx
		jz		end1			// Exit at this point if nothing to do
		add		ecx,ecx
		add		ecx,ecx
		add		ecx,ecx			// ecx now counts bytes
		mov		eax,[s1]		// {{GCC I "a" (s1)}} {{GCC DEL}}
		mov		ebx,[s2]		// {{GCC I "b" (s2)}} {{GCC DEL}}
		add		eax,ecx
		add		ebx,ecx
		neg		ecx				// Now [eax+ecx] etc. point to the 1st element

	top1:
		movdqa	xmm0,[eax+ecx]
		pcmpeqb	xmm0,[ebx+ecx]
		pmovmskb	edx,xmm0	// Unfortunately this has a large latency on the P4 (7 clocks...)
		sub		edx,65535		// Leaves edx zero if all 16 bytes were equal
		jnz		end1
		add		ecx,16
		jnz		top1

	end1:
		mov		[retVal],edx	// {{GCC O "=m" (retVal)}} {{GCC DEL}} {{GCC INS "mov %%edx,%0\n\t"}}
	}

	return retVal;
}

// Generate all combinations of Fitch algorithms using the X-macro technique.
// I *believe* this will ultimately be more maintainable than maintaining separate blobs of code
// for each different combination of switches.
#include "generate_all_fitch_algos_b2_w16_sse2asm.inc"
