// Turn off warnings about unreferenced parameters on known compilers.  To cater to all compilers you need:
//
// #include "nowarnunusedparams.h"
// double some_function(int a_used_arg, NOTUSED(char *a_potentially_unused_arg)) { ... }
// #include "warnpop.h"
//
// for **each** function.  (Borland's "#pragma argsused" only applies to the next function, so this
// can't be used around a block of functions unfortunately.)
//
// We need to use the X-macro technique here since a few compilers (e.g. MSVC, Borland) use #pragmas
// to turn on/off warnings and it's not possible to generate a #pragma line with a preprocessor macro.  :(
#ifdef __BORLANDC__
#pragma argsused
#endif	// __BORLANDC__
#ifdef _MSC_VER
#pragma warning(push)
#pragma warning(disable : 4100)		// "Unreferenced formal parameter"
#endif	// _MSC_VER
#ifdef __GNUC__
#define PARAMNOTUSED(x) x __attribute((unused))
#else	// not __GNUC__
#define PARAMNOTUSED(x) x
#endif	// not __GNUC__
