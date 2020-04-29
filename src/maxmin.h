#ifndef __MAXMIN_H
#define __MAXMIN_H
// An external configuration process should #define INLINE to be either "inline" (if it is supported)
// or blank (in the increasingly unlikely event that it is not).  But if it doesn't:
#ifndef INLINE
// Assume we do have "inline" ability -- that way if we do, everything works happily and quickly; if we
// don't, we get compilation errors which should point out what we need to do.  Would be nice if we
// could warn on this, but there's no portable equivalent of #error for warnings.
#define INLINE inline
#define __MAXMIN_H_WE_DEFINED_INLINE
#endif	// not INLINE

// Believe it or not, MSVC++'s C compiler does not recognise the inline keyword, even though its C++
// compiler does!  Actually it makes no attempt to be C99 compliant at all, and "inline" is a C99-ism.
// But all's not lost, since empirically on either /O1 or /O2 the MSVC++ 2008 will inline these routines anyway.

static INLINE int iMax(int x, int y) {
	return x > y ? x : y;
}

static INLINE int iMin(int x, int y) {
	return x < y ? x : y;
}

static INLINE unsigned uMax(unsigned x, unsigned y) {
	return x > y ? x : y;
}

static INLINE unsigned uMin(unsigned x, unsigned y) {
	return x < y ? x : y;
}

static INLINE long lMax(long x, long y) {
	return x > y ? x : y;
}

static INLINE long lMin(long x, long y) {
	return x < y ? x : y;
}

static INLINE unsigned long ulMax(unsigned long x, unsigned long y) {
	return x > y ? x : y;
}

static INLINE unsigned long ulMin(unsigned long x, unsigned long y) {
	return x < y ? x : y;
}

// Don't bother with separate routines for float, since floats are promoted to doubles
// and can implicitly convert back anyway.
static INLINE double dblMax(double x, double y) {
	return x > y ? x : y;
}

static INLINE double dblMin(double x, double y) {
	return x < y ? x : y;
}

// Clean up after ourselves if necessary...
#ifdef __MAXMIN_H_WE_DEFINED_INLINE
#undef INLINE
#endif	// not __MAXMIN_H_WE_DEFINED_INLINE
#endif	// #include guard
