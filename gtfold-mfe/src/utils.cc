#include "utils.h"
#include <string.h>
#include <string>

char baseToDigit(const char* base) {
	if (!strcmp(base, "A")) {
		return '1';
	}
	if (!strcmp(base, "C")) {
		return '2';
	}
	if (!strcmp(base, "G")) {
		return '3';
	}
	if (!strcmp(base, "U")) {
		return '4';
	}
	if (!strcmp(base, "N")) {
		return '5';
	}
	return (char) '0';
}

unsigned char encode(char base) {
	switch(base) {
		case 'A':
		case 'a':
			return BASE_A;
		case 'C':
		case 'c':
			return BASE_C;
		case 'G':
		case 'g':
			return BASE_G;
		case 'U':
		case 'u':
		case 'T':
		case 't':
			return BASE_U;
		case 'N':
		case 'n':
		case 'R':
		case 'r':
		case 'Y':
		case 'y':
		case 'M':
		case 'm':
		case 'K':
		case 'k':
		case 'S':
		case 's':
		case 'W':
		case 'w':
		case 'B':
		case 'b':
		case 'D':
		case 'd':
		case 'H':
		case 'h':
		case 'V':
		case 'v':
			return BASE_A;
		default:
			return 'X';
	}
}

/**
 * Identical to getBase but returns 'N' instead of 'A' for non-standard bases
 * (a,c,g,u)
 */
unsigned char getBase1(const char* base) {
	if (!strcmp(base, "A") ||
		!strcmp(base, "a"))
		return BASE_A;

	if (!strcmp(base, "C") ||
		!strcmp(base, "c"))
		return BASE_C;

	if (!strcmp(base, "G") ||
		!strcmp(base, "g"))
		return BASE_G;

	if (!strcmp(base, "U") ||
		!strcmp(base, "u") ||
		!strcmp(base, "T") ||
		!strcmp(base, "t"))
		return BASE_U;

	if (!strcmp(base, "N") ||
		!strcmp(base, "n") ||
		!strcmp(base, "R") ||
		!strcmp(base, "r") ||
		!strcmp(base, "Y") ||
		!strcmp(base, "y") ||
		!strcmp(base, "M") ||
		!strcmp(base, "m") ||
		!strcmp(base, "K") ||
		!strcmp(base, "k") ||
		!strcmp(base, "S") ||
		!strcmp(base, "s") ||
		!strcmp(base, "W") ||
		!strcmp(base, "w") ||
		!strcmp(base, "B") ||
		!strcmp(base, "b") ||
		!strcmp(base, "D") ||
		!strcmp(base, "d") ||
		!strcmp(base, "H") ||
		!strcmp(base, "h") ||
		!strcmp(base, "V") ||
		!strcmp(base, "v"))
		return 'N';

	return 'X';
}
