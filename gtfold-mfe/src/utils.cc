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

unsigned char getBase(const char* base) {
	//cout << base;
	if (!strcmp(base, "A") || !strcmp(base, "a")) {
		//cout << "1";
		return BASE_A;
	}
	if (!strcmp(base, "C") || !strcmp(base, "c")) {
		//cout << "2";
		return BASE_C;
	}
	if (!strcmp(base, "G") || !strcmp(base, "g")) {
		//cout << "3";
		return BASE_G;
	}
	if (!strcmp(base, "U") || !strcmp(base, "u") || !strcmp(
			base, "T") || !strcmp(base, "t")) {
		//cout << "4";
		return BASE_U;
	}
	if (!strcmp(base, "N") || !strcmp(base, "n")||!strcmp(base, "R")|| !strcmp(base, "r")|| !strcmp(base, "Y")|| !strcmp(base, "y")|| !strcmp(base, "M")|| !strcmp(base, "m")|| !strcmp(base, "K")|| !strcmp(base, "k")|| !strcmp(base, "S")|| !strcmp(base, "s")|| !strcmp(base, "W")|| !strcmp(base, "w")|| !strcmp(base, "B")|| !strcmp(base, "b")|| !strcmp(base, "D")|| !strcmp(base, "d")|| !strcmp(base, "H")|| !strcmp(base, "h")|| !strcmp(base, "V")|| !strcmp(base, "v")) {
		//cout << "4";
		return BASE_A;
	}
	return 'X';
}

unsigned char getBase1(const char* base) {
	//cout << base;
	if (!strcmp(base, "A") || !strcmp(base, "a")) {
		//cout << "1";
		return BASE_A;
	}
	if (!strcmp(base, "C") || !strcmp(base, "c")) {
		//cout << "2";
		return BASE_C;
	}
	if (!strcmp(base, "G") || !strcmp(base, "g")) {
		//cout << "3";
		return BASE_G;
	}
	if (!strcmp(base, "U") || !strcmp(base, "u") || !strcmp(
																			base, "T") || !strcmp(base, "t")) {
		//cout << "4";
		return BASE_U;
	}
	if (!strcmp(base, "N") || !strcmp(base, "n")|| !strcmp(base, "R")|| !strcmp(base, "r")|| !strcmp(base, "Y")|| !strcmp(base, "y")|| !strcmp(base, "M")|| !strcmp(base, "m")|| !strcmp(base, "K")|| !strcmp(base, "k")|| !strcmp(base, "S")|| !strcmp(base, "s")|| !strcmp(base, "W")|| !strcmp(base, "w")|| !strcmp(base, "B")|| !strcmp(base, "b")|| !strcmp(base, "D")|| !strcmp(base, "d")|| !strcmp(base, "H")|| !strcmp(base, "h")|| !strcmp(base, "V")|| !strcmp(base, "v")) {
		//cout << "4";
		return 'N';
	}
	return 'X';
}
