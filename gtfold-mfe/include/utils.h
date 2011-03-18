#ifndef _UTILS_H_
#define _UTILS_H_

#include "constants.h"


#define MIN(X,Y) ((X) < (Y) ? (X) : (Y))
#define MAX(X,Y) ((X) > (Y) ? (X) : (Y))

char baseToDigit(const char* base) ;
unsigned char encode(char base);
unsigned char getBase1(const char* base) ;

#endif
