/*
 GTfold: compute minimum free energy of RNA secondary structure
 Copyright (C) 2008  David A. Bader
 http://www.cc.gatech.edu/~bader

 This program is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef _ALGORITHMS_H
#define _ALGORITHMS_H

#include "main-c.h"

#define MIN(X,Y) ((X) < (Y) ? (X) : (Y))
#define MAX(X,Y) ((X) > (Y) ? (X) : (Y))

/* It returns zero if (i,j) can not make a base pair otherwise return 1. The possible base pairs are (A,U), (U,A), (C,G), (G,C), (G,U) and (U,G). Note that the following condition results into 1 for the allowed base pairs .*/
#define checkPair(i, j) (((((i)-(j)) % 2) == 1 || (((i)-(j)) % 2)== -1) && (!( ((i)==BASE_A && (j)==BASE_C) || ((i)==BASE_C && (j)==BASE_A) )))
/* Non GC penalty or AU penalty. Returns a constant penalty if the base pair is not GC or CG */
#define auPen(i, j) ((( (i)==BASE_U || (j)==BASE_U ) && ( (i)==BASE_A || (i)==BASE_G || (j)==BASE_A || (j)==BASE_G )) ? auend : 0)

extern int *constraints;

#ifdef __cplusplus
extern "C" {
#endif
	void initTables(int len);
	
	void init_chPair();
	int update_chPair(int i, int j);
#ifdef __cplusplus
}
#endif

//void progArrays();

#ifdef __cplusplus
extern "C" {
#endif
	int calculate(int len, int **forceList, int **prohibitList, int forcelen,
			int prohibitlen);
	enum BOOL {
		FALSE, TRUE
	};
	extern enum BOOL ILSA;
	extern enum BOOL NOISOLATE;
#ifdef __cplusplus
}
#endif

int checkSS(int i, int j);

void calcVBI(int i, int j);
void calcVM(int i, int j);
void calcWM(int i, int j);
void calcW(int j);

void calcVWM(int i, int j, int vbiij, int vmij);
void calcVBIVMVWM(int i, int j);
void calcVBIVMVWM2(int i, int j);

void calcVBIS(int i, int j);
void calcVMVWM(int i, int j);
int eS(int i, int j);
int eH(int i, int j);
int eL(int i, int j, int ip, int jp);

#endif

