#include <stdio.h>
#include <stdlib.h>

//#define PASTE(x, y) PASTE2(x, y)
//#define PASTE2(x, y) x ## y
#define REPORT_SIZEOF(t) printf("sizeof (" #t ")=%d\n", sizeof (t));

int main(int argc, char **argv) {
	REPORT_SIZEOF(char);
	REPORT_SIZEOF(short);
	REPORT_SIZEOF(int);
	REPORT_SIZEOF(long);
	REPORT_SIZEOF(long long);
	REPORT_SIZEOF(void *);
	return 0;
}
