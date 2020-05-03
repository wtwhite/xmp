#ifndef __ARENA_H
#define __ARENA_H
// A simple arena data structure for performing LIFO memory allocations quickly.
// ArenaInit() calls ALIGNED_MALLOC(), so provided that all allocations are multiples
// of BYTESPERBLOCK, all returned memory blocks are aligned.

#include <stdlib.h>			// For size_t
#include <assert.h>
#include <stdint.h>			// Needed for uintptr_t
#include <string.h>
#include "switches.h"		// Needed for BYTESPERBLOCK
#include "alignedalloc.h"	// Needed for ALIGNED_MALLOC(), ALIGNED_FREE()
#include <dbgprint.h>

struct arena {
	unsigned char *cur;
	unsigned char *base;
	size_t size;
};

// An external configuration process should #define INLINE to be either "inline" (if it is supported)
// or blank (in the increasingly unlikely event that it is not).

static INLINE void ArenaCreate(struct arena *a, size_t size) {
	DBGPRINT3("ArenaCreate(a=%p, size=%u) called.\n", a, size);
	assert(size % BYTESPERBLOCK == 0);
	ALIGNED_MALLOC(a->base, size);
	assert((uintptr_t) a->base % BYTESPERBLOCK == 0);
	a->cur = a->base;
	a->size = size;
}

static INLINE void ArenaDestroy(struct arena *a) {
	DBGPRINT2("ArenaDestroy(a=%p) called.\n", a);
	ALIGNED_FREE(a->base);
}

static INLINE void *ArenaAllocate(struct arena *a, size_t size) {
	unsigned char *p = a->cur;
	DBGPRINT3("ArenaAllocate(a=%p, size=%u) called.\n", a, size);
	assert(size % BYTESPERBLOCK == 0);
	assert(a->cur + size <= a->base + a->size);
	a->cur += size;

#ifdef DEBUG
	// Mimic MS's debug allocation strategy
	memset(p, 0xCD, size);
#endif	// DEBUG

	return p;
}

static INLINE void *ArenaGetPos(struct arena *a) {
	DBGPRINT3("ArenaGetPos(a=%p) called, returning %p.\n", a, a->cur);
	return a->cur;
}

static INLINE void ArenaSetPos(struct arena *a, void *p) {
	DBGPRINT3("ArenaSetPos(a=%p, p=%p) called.\n", a, p);
	assert((unsigned char *) p >= a->base && (unsigned char *) p <= a->base + a->size);
	assert((uintptr_t) p % BYTESPERBLOCK == 0);
	a->cur = p;
}
#endif	// #include guard
