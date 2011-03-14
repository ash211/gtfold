#ifndef _GLOBAL_H_
#define _GLOBAL_H_

#include "constants.h"

extern unsigned char *RNA1; 
extern unsigned char *RNA; 
extern int *structure; 
extern int VERBOSE;

extern unsigned int chPairKey;

// The possible base pairs are (A,U), (U,A), (C,G), (G,C), (G,U) and (U,G). 
#define checkPair(i, j) (((((i)-(j)) % 2) == 1 || (((i)-(j)) % 2)== -1) && (!( ((i)==BASE_A && (j)==BASE_C) || ((i)==BASE_C && (j)==BASE_A) )))

inline int allowedPairs(int a, int b) ; 

#ifdef __cplusplus
extern "C" {
#endif
void init_global_params(int len);
void free_global_params();
void printSequence(int len); 
void printStructure(int len); 
#ifdef __cplusplus
}
#endif

void init_checkPair(); 
int  update_checkPair(int i, int j);

#endif
