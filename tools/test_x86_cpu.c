#include <stdio.h>
#include <assert.h>

int verbose = 0;
int quiet = 0;
static int test_result;

#ifdef _MSC_VER
#include <windows.h>

#define PRE __try {
#define POST														\
	}																\
	__except (EXCEPTION_EXECUTE_HANDLER) {							\
		if (verbose) printf("Exception 0x%x occurred!\n", GetExceptionCode());	\
		test_result = 0;											\
	}
#endif	// _MSC_VER

#ifdef __unix__
#include <signal.h>
#include <setjmp.h>

#define PRE															\
	{																\
		static struct sigaction sa, oldsa;							\
		sa.sa_handler = &sig_handler;								\
		sigemptyset(&sa.sa_mask);									\
		/* SA_NODEFER is important -- otherwise SIGILL remains blocked (because the handler calls longjmp() instead of */		\
		/* returning normally), causing program abort on the **2nd** invalid opcode. */		\
		/* Could equivalently use sigsetjmp()+siglongjmp(), or sigprocmask() -- all require POSIX.1-2001. */						\
		sa.sa_flags = SA_NODEFER | SA_RESETHAND;					\
		sigaction(SIGILL, &sa, &oldsa);								\
		if (!setjmp(jbuf)) {
#define POST														\
		}															\
	}

jmp_buf jbuf;
void sig_handler(int sig) {
	if (verbose) printf("Signal %d occurred!\n", sig);
	test_result = 0;
	longjmp(jbuf, 1);
}
#endif	// __unix__

//void mmx_test() {
//	__asm {
//		pand	mm0,mm0
//		emms
//	};
//}
//
//void sse_test() {
//	__asm {
//		movaps	xmm0,xmm0
//	};
//}
//
//void sse2_test() {
//	__asm {
//		pand	xmm0,xmm0
//	};
//}

struct entry {
	char *name;
	char *code;
} tests[] = {
	{ "nop",   "\x90\xC3" },					// nop; ret (negative control)
	{ "crash", "\x0F\x0B\xC3" },				// 0F 0B is one of 2 recommended byte sequences to generate "invalid instruction" exception on x86 (positive control)
	{ "mmx",   "\x0F\xDB\xC0\x0F\x77\xC3" },	// pand mm0,mm0; emms; ret
	{ "sse",   "\x0F\x28\xC0\xC3" },			// movaps xmm0,xmm0; ret
	{ "sse2",  "\x66\x0F\xDB\xC0\xC3" }			// pand xmm0,xmm0; ret
};

int test(char *name) {
	int i;
	for (i = 0; i < sizeof tests / sizeof tests[0]; ++i) {
		if (!strcmp(tests[i].name, name)) {
			if (verbose) printf("About to test '%s'...\n", name);
			test_result = 1;
			
			PRE
			((void (*)()) tests[i].code)();
			POST
			
			if (verbose) printf("Finished test.\n");
			if (!quiet) {
				printf("%s %s\n", name, test_result ? "PASS" : "FAIL");
			}
			return test_result;
		}
	}
	
	printf("Could not find test '%s'.\n", name);
	return -1;
}

int main(int argc, char **argv) {
	int i;
	int first_failure = 0;
	
	// Only 2 possible options and they're mutually exclusive.
	if (!strcmp(argv[1], "-v") || !strcmp(argv[1], "--verbose")) {
		verbose = 1;
		++argv;
		--argc;
	} else if (!strcmp(argv[1], "-q") || !strcmp(argv[1], "--quiet")) {
		quiet = 1;
		++argv;
		--argc;
	}
	
	for (i = 1; i < argc; ++i) {
		if (!test(argv[i]) && !first_failure) {
			first_failure = i;
		}
	}
	
	if (verbose) printf("Finished all tests.  First failure on test #%d.\n", first_failure);
	return first_failure;
}
