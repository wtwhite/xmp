PROJECT(fastdnamp)

CMAKE_MINIMUM_REQUIRED(VERSION 2.6)

INCLUDE(CheckCSourceCompiles)		# Used by several things

IF(NOT DEFINED USEFUL_LIBS_INCLUDE_DIR)
	IF(IS_DIRECTORY ${CMAKE_SOURCE_DIR}/../useful_libs)
		SET(USEFUL_LIBS_INCLUDE_DIR ${CMAKE_SOURCE_DIR}/../useful_libs/include CACHE FILEPATH "Include directory for useful_libs")
	ENDIF(IS_DIRECTORY ${CMAKE_SOURCE_DIR}/../useful_libs)
ENDIF(NOT DEFINED USEFUL_LIBS_INCLUDE_DIR)

IF(NOT DEFINED USEFUL_LIBS_INCLUDE_DIR)
	IF(IS_DIRECTORY ${CMAKE_SOURCE_DIR}/../../useful_libs)
		SET(USEFUL_LIBS_INCLUDE_DIR ${CMAKE_SOURCE_DIR}/../../useful_libs/include CACHE FILEPATH "Include directory for useful_libs")
	ENDIF(IS_DIRECTORY ${CMAKE_SOURCE_DIR}/../../useful_libs)
ENDIF(NOT DEFINED USEFUL_LIBS_INCLUDE_DIR)

#HACK: Do the same thing for libs?  Harder though.

#HACK: Need to incorporate paths to my libraries better somehow...
MESSAGE(CMAKE_SYSTEM_NAME=${CMAKE_SYSTEM_NAME})
IF(WIN32)
	SET(OS win)
	IF(NOT DEFINED USEFUL_LIBS_INCLUDE_DIR)
		SET(USEFUL_LIBS_INCLUDE_DIR H:/cprogs/useful_libs/include CACHE FILEPATH "Include directory for useful_libs")
	ENDIF(NOT DEFINED USEFUL_LIBS_INCLUDE_DIR)
	
	IF(NOT DEFINED USEFUL_LIBS_LIB_DIR)
		SET(USEFUL_LIBS_LIB_DIR H:/cprogs/useful_libs/lib/prebuilt/win32/msvc9-dllcrt-release CACHE FILEPATH "Library directory for useful_libs")
	ENDIF(NOT DEFINED USEFUL_LIBS_LIB_DIR)
ELSE(WIN32)
	IF(NOT DEFINED USEFUL_LIBS_INCLUDE_DIR)
		SET(USEFUL_LIBS_INCLUDE_DIR ~/H/cprogs/useful_libs/include CACHE FILEPATH "Include directory for useful_libs")
	ENDIF(NOT DEFINED USEFUL_LIBS_INCLUDE_DIR)
	
	STRING(COMPARE EQUAL ${CMAKE_SYSTEM_NAME} Linux LINUX)
	IF(LINUX)
		SET(OS linux)
		
		#HACK: Empirically "x86" on Windows, "i686" on Linux...
		STRING(REGEX MATCH "^[ix].*86\$" INTEL ${CMAKE_SYSTEM_PROCESSOR})
		IF(INTEL)
			IF(NOT DEFINED USEFUL_LIBS_LIB_DIR)
				SET(USEFUL_LIBS_LIB_DIR ~/H/cprogs/useful_libs/lib/prebuilt/linux-x86 CACHE FILEPATH "Library directory for useful_libs")
			ENDIF(NOT DEFINED USEFUL_LIBS_LIB_DIR)
			MESSAGE("Set link directories for Linux x86.")
		ENDIF(INTEL)
	ELSE(LINUX)
		# Needed for timer_${OS}.c
		SET(OS dummy)
	ENDIF(LINUX)
ENDIF(WIN32)

IF(NOT DEFINED USEFUL_LIBS_LIB_DIR)
	MESSAGE(SEND_ERROR "You must specify USEFUL_LIBS_LIB_DIR.")
ENDIF(NOT DEFINED USEFUL_LIBS_LIB_DIR)
#SET(USEFUL_LIBS_LIB_DIR ${USEFUL_LIBS_LIB_DIR} CACHE FILEPATH "Library directory for useful_libs")		#HACK: Hopefully this puts it in the cache...

INCLUDE_DIRECTORIES(${USEFUL_LIBS_INCLUDE_DIR})
LINK_DIRECTORIES(${USEFUL_LIBS_LIB_DIR})

# Generated C files get added to the build dir, so we need to explicitly add
# the source dir to the include path.
INCLUDE_DIRECTORIES(${INCLUDE_DIRECTORIES} ${CMAKE_SOURCE_DIR})

# Big-endian or little-endian?
INCLUDE(TestBigEndian)
TEST_BIG_ENDIAN(BIGENDIAN)
#MESSAGE(BIGENDIAN=${BIGENDIAN})
IF(BIGENDIAN)
	ADD_DEFINITIONS(-DBIGENDIAN)
ELSE(BIGENDIAN)
	ADD_DEFINITIONS(-DLITTLEENDIAN)
ENDIF(BIGENDIAN)

# We need the stdint.h header from C99 for types like uint32_t, but some environments (e.g. MSVC++) don't include it.
# If it isn't on this system, use a portable one from http://www.azillionmonkeys.com/qed/pstdint.h.  Because we
# include the build dir in the include path, we can just use "#include <stdint.h>" as usual in the source code.
# Note to self: Lucky I chanced upon CheckIncludeFiles, since cmake's documentation of FIND_FILE() is bogus.
# None of the mentioned variables (e.g. CMAKE_SYSTEM_INCLUDE_PATH) contains the path
# "C:/Program Files/Microsoft Visual Studio 9.0/VC", or that with "include" at the end, even though mysteriously
# FIND_FILE(STDIO_H stdio.h) DOES find the stdio.h file there.  Ugh.  (I would have liked to use FIND_FILE() directly,
# but couldn't figure out how to make it search *only* the system include paths.)
INCLUDE(CheckIncludeFiles)
CHECK_INCLUDE_FILES(stdint.h HAS_STDINT_H)
IF(NOT HAS_STDINT_H)
	EXECUTE_PROCESS(COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_SOURCE_DIR}/pstdint.h ${CMAKE_BINARY_DIR}/stdint.h)
ENDIF(NOT HAS_STDINT_H)

CHECK_C_SOURCE_COMPILES("inline int f() { return 42; } int main(int argc, char **argv) { return 0; }" HAS_INLINE_KEYWORD)
IF(HAS_INLINE_KEYWORD)
	SET(INLINE_KEYWORD "inline" CACHE STRING "'inline' if your C compiler supports the 'inline' keyword, otherwise blank.")
ELSE(HAS_INLINE_KEYWORD)
	SET(INLINE_KEYWORD "" CACHE STRING "'inline' if your C compiler supports the 'inline' keyword, otherwise blank.")
ENDIF(HAS_INLINE_KEYWORD)

CHECK_C_SOURCE_COMPILES("void __fastcall f() {} int main(int argc, char **argv) { return 0; }" HAS_FASTCALL_KEYWORD)
IF(HAS_FASTCALL_KEYWORD)
	SET(FASTCALL_KEYWORD "__fastcall" CACHE STRING "'__fastcall' if your C compiler supports the '__fastcall' keyword, otherwise blank.")
ELSE(HAS_FASTCALL_KEYWORD)
	CHECK_C_SOURCE_COMPILES("void __attribute__((regparm(3))) f() {} int main(int argc, char **argv) { return 0; }" HAS_REGPARM_ATTRIBUTE)
	IF(HAS_REGPARM_ATTRIBUTE)
		SET(FASTCALL_KEYWORD "__attribute__((regparm(3)))" CACHE STRING "'__fastcall' if your C compiler supports the '__fastcall' keyword, otherwise blank.")
	ELSE(HAS_REGPARM_ATTRIBUTE)
		SET(FASTCALL_KEYWORD "" CACHE STRING "'__fastcall' if your C compiler supports the '__fastcall' keyword, otherwise blank.")
	ENDIF(HAS_REGPARM_ATTRIBUTE)
ENDIF(HAS_FASTCALL_KEYWORD)

SET(SRCS
	fastdnamp.c
	common.c
	bound.c
	distbasesbound.c
	incompatbound.c
	mstbound.c
	hittingsetbound.c
	yanbader.c
	seq_b1.c
	seq_b2_w1_basic.c
	seq_b2_w4_fastc.c
	sitepaircosts.c
	ordertaxa.c
	partbound.c
	maketempfile.c
	greedytree.c
	tbr.c
)

# Create parallel MPI version?
OPTION(TARGETMULTI "Produce a parallel MPI version of the software")
IF(TARGETMULTI)
	# We can build MPI programs one of two ways -- either use a "plain" compiler and get the FindMPI module to
	# look for the MPI headers and libraries on the local system, or rely on the user supplying a "special" compiler
	# via CMAKE_C_COMPILER that effectively wraps the compiler (e.g. this might be a script called mpicc or similar).
	# The latter approach is more likely to work on "tricky" setups since FindMPI's heuristics are pretty weak; but
	# the downside is that various other tests performed by CMake will also use this compiler, so it must be
	# able to compile regular non-MPI code that can be run without mpirun/mpiexec.
	# Rather than ask the user to explicitly specify which approach they want, we can just assume the latter and
	# test-compile a small program that calls MPI functions -- if it fails, assume the former and call FindMPI.
	CHECK_C_SOURCE_COMPILES("
#include \"mpi.h\"
int main(int argc, char **argv) { MPI_Init(&argc, &argv); MPI_Finalize(); return 0; }
"
		MPI_CONFIGURED_IMPLICITLY
	)
	IF(NOT MPI_CONFIGURED_IMPLICITLY)
		INCLUDE(FindMPI)
		IF(NOT MPI_FOUND)
			MESSAGE(SEND_ERROR "Cannot find MPI libraries.  Please specify MPI_INCLUDE_PATH and MPI_LIBRARY manually.")
		ENDIF(NOT MPI_FOUND)
		INCLUDE_DIRECTORIES(${MPI_INCLUDE_PATH})
	ENDIF(NOT MPI_CONFIGURED_IMPLICITLY)
	SET(SRCS ${SRCS} mpiboss.c mpiworker.c)
	
	# Use MPI synchronous sends (MPI_Issend()), or the unproven but probably faster standard-mode sends?
	OPTION(SYNCHRONOUSSENDS "Use safe but slow synchronous sends (MPI_Issend()) for the parallel MPI version")
ENDIF(TARGETMULTI)

# Add in the b2_w8 variation on 64-bit machines.
INCLUDE(CheckTypeSize)
#HACK: You may need to enable one of the following switches to get an efficient binary on a 64-bit computer...
#SET(CMAKE_REQUIRED_FLAGS -m64)		#HACK: Enable 64-bit binaries with gcc
#SET(CMAKE_REQUIRED_FLAGS -q64)		#HACK: Enable 64-bit binaries with xlc
CHECK_TYPE_SIZE("unsigned long" WORDSIZE BUILTIN_TYPES_ONLY)
IF(WORDSIZE EQUAL 8 OR WORDSIZE GREATER 8)
	MESSAGE("64-bit machine detected.")
	SET(FASTC64FITCH ON CACHE BOOL "Use optimised 64-bit C routines for calculation of Fitch algorithm (recommended for 64-bit architectures, but will probably work (slowly) on 32-bit architectures)")
ELSE(WORDSIZE EQUAL 8 OR WORDSIZE GREATER 8)
	MESSAGE("32-bit machine detected.")
ENDIF(WORDSIZE EQUAL 8 OR WORDSIZE GREATER 8)

# Select FASTCFITCH as a sensible default if nothing else is selected.
IF(NOT BASICFITCH AND NOT FASTCFITCH AND NOT FASTC64FITCH AND NOT SSE2ASMFITCH)
	MESSAGE("Defaulting to FASTCFITCH.")
	SET(FASTCFITCH ON CACHE BOOL "Use optimised C routines for calculation of Fitch algorithm (all architectures)")
ENDIF(NOT BASICFITCH AND NOT FASTCFITCH AND NOT FASTC64FITCH AND NOT SSE2ASMFITCH)

#TODO: If the CPU looks like an x86, probe it to see whether it will run SSE2 code.
OPTION(DEBUG "Produce extra debugging output" ON)
OPTION(BASICFITCH "Use slow C routines for calculation of Fitch algorithm (all architectures)")
OPTION(FASTCFITCH "Use optimised C routines for calculation of Fitch algorithm (all architectures)")
OPTION(FASTC64FITCH "Use optimised 64-bit C routines for calculation of Fitch algorithm (recommended for 64-bit architectures, but will probably work (slowly) on 32-bit architectures)")
OPTION(SSE2ASMFITCH "Use SSE2 assembly routines for calculation of Fitch algorithm (Requires Pentium 4 on Windows)")
SET(SITESPERBYTE 2 CACHE STRING "Number of sites to pack into a single byte (1 or 2)")

# I guess this won't compile on Linux...
IF(WIN32)
	MACRO(ENSURE_ASM_EXISTS BASENAME)
		# Do nothing.  Windows versions of C files containing inline assembly are the base files.
	ENDMACRO(ENSURE_ASM_EXISTS BASENAME)
	ADD_DEFINITIONS(-D_CRT_SECURE_NO_WARNINGS)	# Get rid of annoying MSVC++ warnings about fopen() et al. being "deprecated"(!)
	SET(TIMER_C timer_win.c)
	SET(ELAPSEDTIME_C elapsedtime_win.c)
ENDIF(WIN32)

# Need to generate the gas-syntax inline-assembly files from the Intel-syntax ones.
IF(UNIX)
	MACRO(CONVERT_ASM FROM TO)
		ADD_CUSTOM_COMMAND(
			OUTPUT ${TO}
			COMMAND ${CMAKE_SOURCE_DIR}/../tools/convert_intel_asm_to_gcc.pl ARGS ${CMAKE_CURRENT_SOURCE_DIR}/${FROM} > ${TO}
			DEPENDS
				${CMAKE_CURRENT_SOURCE_DIR}/${FROM}
				${CMAKE_SOURCE_DIR}/../tools/convert_intel_asm_to_gcc.pl
			COMMENT "Generating gas-syntax assembly from ${FROM} for ${TO}..."
		)
	ENDMACRO(CONVERT_ASM FROM TO)
	
	MACRO(ENSURE_ASM_EXISTS BASENAME)
		CONVERT_ASM(${BASENAME}_win.c ${BASENAME}_${OS}.c)
		
		# Seems recent versions of gcc (and possibly old versions too) require -march to be set
		# to recognise MMX register names like %mm1, %mm2 etc.
		STRING(REPLACE -march=pentium4 "" CMAKE_C_FLAGS "${CMAKE_C_FLAGS}")
		SET(CMAKE_C_FLAGS ${CMAKE_CFLAGS} -march=pentium4)
	ENDMACRO(ENSURE_ASM_EXISTS BASENAME)
	
	# Do we support the usual functions for creating timers?  (HACK: Should check for sigaction() too...)
	CHECK_C_SOURCE_COMPILES("
#include <unistd.h>
int main(int argc, char **argv) { alarm(42); return 0; }
"
		HAS_POSIX_TIMER
	)
	IF(HAS_POSIX_TIMER)
		SET(TIMER_C timer_posix.c)
	ENDIF(HAS_POSIX_TIMER)
	
	CHECK_C_SOURCE_COMPILES("
#include <time.h>
#include <sys/time.h>
int main(int argc, char **argv) { gettimeofday(NULL, NULL); return 0; }
"
		HAS_GETTIMEOFDAY
	)
	IF(HAS_GETTIMEOFDAY)
		SET(ELAPSEDTIME_C elapsedtime_posix.c)
	ENDIF(HAS_GETTIMEOFDAY)
ENDIF(UNIX)

IF(NOT DEFINED TIMER_C)
	SET(TIMER_C timer_dummy.c)		# These routines do nothing.  We can live without periodic updates.
ENDIF(NOT DEFINED TIMER_C)
IF(NOT DEFINED ELAPSEDTIME_C)
	MESSAGE(SEND_ERROR "Could not find a way to measure time durations on your platform!")
ENDIF(NOT DEFINED ELAPSEDTIME_C)

SET(SRCS ${SRCS} ${TIMER_C} ${ELAPSEDTIME_C})

SET(NUM_IMPLS 0)
IF(SSE2ASMFITCH)
	SET(BYTESPERBLOCK 16 CACHE STRING "Number of bytes per block" FORCE)
	#HACK: The following is not yet compiling properly on Linux.
	SET(SRCS ${SRCS} seq_b2_sse2asm_${OS}.c)
	ENSURE_ASM_EXISTS(seq_b2_sse2asm)
	MATH(EXPR NUM_IMPLS ${NUM_IMPLS}+1)
ENDIF(SSE2ASMFITCH)
IF(FASTCFITCH)
	SET(BYTESPERBLOCK 4 CACHE STRING "Number of bytes per block")
	MATH(EXPR NUM_IMPLS ${NUM_IMPLS}+1)
ENDIF(FASTCFITCH)
IF(FASTC64FITCH)
	SET(BYTESPERBLOCK 8 CACHE STRING "Number of bytes per block")
	SET(SRCS ${SRCS} seq_b2_w8_fastc64.c)
	MATH(EXPR NUM_IMPLS ${NUM_IMPLS}+1)
ENDIF(FASTC64FITCH)
IF(BASICFITCH)
	MATH(EXPR NUM_IMPLS ${NUM_IMPLS}+1)
ENDIF(BASICFITCH)
IF(NOT DEFINED BYTESPERBLOCK)
	# This could happen e.g. for BASICFITCH, for which no default block size is given.
	MESSAGE(SEND_ERROR "You must specify a value for BYTESPERBLOCK.")
ENDIF(NOT DEFINED BYTESPERBLOCK)

OPTION(STOPEARLY "Stop computing the Fitch algorithm as soon as the score exceeds the bound.  Usually gives a small speedup." TRUE)
OPTION(FASTWEIGHT1FITCH "Optimise Fitch calculations for weight-1 columns.  Usually gives a moderate speedup." TRUE)
OPTION(WEIGHT1MULTIPLY "Use a multiply operation to sum weight-1 column scores for the Fitch algorithm.  On many architectures this is faster than repeated additions." TRUE)
#HACK: MEASURE should default to FALSE in a production environment.
OPTION(MEASURE "Calculate and output a variety of diagnostic measurements and statistics at the end of analysis." TRUE)

# Sanity check
IF(NOT NUM_IMPLS EQUAL 1)
	MESSAGE(SEND_ERROR "You must specify exactly one of SSE2ASMFITCH, FASTCFITCH, FASTC64FITCH and BASICFITCH.")
ENDIF(NOT NUM_IMPLS EQUAL 1)

IF(MEASURE)
	SET(SRCS ${SRCS} measure.c)
ENDIF(MEASURE)

#DEBUG
MESSAGE(CMAKE_SYSTEM_PROCESSOR=${CMAKE_SYSTEM_PROCESSOR})

# Produce switches.h in the build dir, which contains #defines for configured options like DEBUG and SSE2ASMFITCH.
# It is generated at configure time, so it will be picked up properly by build-time dependency checking.
CONFIGURE_FILE(${CMAKE_SOURCE_DIR}/switches.h.in ${CMAKE_BINARY_DIR}/switches.h)
INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR})

ADD_EXECUTABLE(fastdnamp ${SRCS})
TARGET_LINK_LIBRARIES(fastdnamp useful)

IF(TARGETMULTI AND NOT MPI_CONFIGURED_IMPLICITLY)
	TARGET_LINK_LIBRARIES(fastdnamp ${MPI_LIBRARIES})
ENDIF(TARGETMULTI AND NOT MPI_CONFIGURED_IMPLICITLY)

OPTION(MPITRACE "Turn on MPE tracing for MPI.  Only useful when TARGETMULTI is turned on.  Requires MPE libraries (e.g these come with MPICH).")
IF(MPITRACE)
	TARGET_LINK_LIBRARIES(fastdnamp lmpe mpe)
ENDIF(MPITRACE)

# MakeTempFilename() calls gethostname() which requires Winsock on Windows.
# Also Windows' GetTimeNow() calls timeGetTime(), which requires winmm.lib.
IF(WIN32)
	TARGET_LINK_LIBRARIES(fastdnamp wsock32 winmm)
ENDIF(WIN32)
