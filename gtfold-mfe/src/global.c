#include <stdio.h>
#include <stdlib.h>

#include "global.h"

unsigned char *RNA1; 
unsigned char *RNA; 
int *structure; 
unsigned int chPairKey;


void init_global_params(int len)
{
	RNA = (unsigned char *) malloc((len+1)* sizeof(unsigned char));
	if (RNA == NULL) {
		perror("Cannot allocate variable 'RNA'");
		exit(-1);
	}
	RNA1 = (unsigned char *) malloc((len+1) * sizeof(unsigned char));
	if (RNA1 == NULL) {
		perror("Cannot allocate variable 'RNA'");
		exit(-1);
	}
	structure = (int *) malloc((len+1) * sizeof(int));
	if (structure == NULL) {
		perror("Cannot allocate variable 'structure'");
		exit(-1);
	}
	
	init_checkPair();
}

void free_global_params()
{
	free(structure);
	free(RNA);
	free(RNA1);
}

void printSequence(int len) 
{
	int i;
	for (i = 1; i <= len; i++) 
	{
		if (RNA1[i] == 0)
			printf("A");
		else if (RNA1[i] == 1)
			printf("C");
		else if (RNA1[i] == 2)
			printf("G");
		else if (RNA1 [i] == 3)
			printf("T");
		else
			printf("N");
	}
	printf("\n");
}

void printStructure(int len) 
{
	int i = 1;
	for (i = 1; i <= len; i++) 
	{
		if (structure[i] > 0 && structure[i] > i)
			printf("(");
		else if (structure[i] > 0 && structure[i] < i)
			printf(")");
		else
			printf(".");
	}
	printf("\n");
}

int canPair(int a, int b)
{
	return (chPairKey & (1 << (((a)<<2) + (b)))); 
}  

void init_checkPair() 
{
	int i, j;
	chPairKey = 0;
	for (i = 0; i < 4; i++)
		for (j = 0; j < 4; j++)
			chPairKey += checkPair(i, j) << ((i << 2) + j);
}

int update_checkPair(int i, int j) 
{
	int r = 0;
	if (!((i >= 0 && i <=3 )&&(j >=0 && j <=3)))
		return r;
	if (!(chPairKey & (1 << ((i << 2) + j))))
	{
		chPairKey += 1 << ((i << 2) + j);	
		r = 1;
	}
	return r;
}
