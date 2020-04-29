#ifndef __ALIGNEDALLOC_H
#define __ALIGNEDALLOC_H

#include "switches.h"

// Simple macros for allocating and freeing memory with a given alignment.
// Needed because Windows and other platforms have different functions for doing this.

#ifdef _MSC_VER
#include <malloc.h>
#define ALIGNED_MALLOC(p, size) (p = _aligned_malloc((size), BYTESPERBLOCK))
#define ALIGNED_FREE(p) _aligned_free((p))
#else	// not _MSC_VER
#include <stdlib.h>
#define ALIGNED_MALLOC(p, size) do { if (posix_memalign(&(p), BYTESPERBLOCK, (size))) assert(!"ALIGNED_MALLOC() failed!"); } while (0)
#define ALIGNED_FREE(p) free((p))
#endif	// not _MSC_VER
#endif	// #include guard
