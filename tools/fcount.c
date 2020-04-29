#include <stdlib.h>
#include <stdio.h>

char *CharName(char c);

int freq[256];

int main(int argc, char **argv) {
	int ch;
	int i;

	while ((ch = getchar()) != EOF) {
		++freq[ch];
	}

	for (i = 0; i < 256; ++i) {
		if (freq[i]) {
			printf("%s: %d\n", CharName((char) i), freq[i]);
		}
	}
}

char *CharName(char c) {
	static char s[2] = { 0, 0 };

	switch (c) {
		case '\n' : return "\\n";
		case '\t' : return "\\t";
		default: s[0] = c; return s;
	}
}
