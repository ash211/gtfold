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

/* AUTHORED by Amrita Mathuriya August 2007 - January 2009. Implemented multiloop energy function, internal loop energy function using heuristic and internal loop speedup algorithm, parallelization, corrected numerous bugs and commented whole GTfold.

 * Amrita: Please note that, in this file same recursion formulas is being calculated in more than one functions.
 * This is done for performance improvement to reduce the redundant computations for various cases. The duplicate codes are not documented again at some places.
 * The arrays to be calculated are VBI, VM, V and WM for every  point (i,j) where j > i. Then W(j) needs to be calculated for j= 1 to N.
 * NOTE that the WM(i,j) can be calculated only after V(i,j) and array VBI and VM should be calculated before V array for point (i,j). So, the order of computation has been kept as VBI, VM, V, WM for any point (i,j)
 * Also Note that, a valid base pair has j > i. Therefore, the portion of the 2D arrays containing j < i is not useful.
 * Minimum size of a hairpin loop is assumed as 3. This assumption is taken into effect at many places.
 * I am not sure of what these eparam values are at various places except for multiloops.
 * */

/* Modified by Sainath Mallidi August 2009 -  "*/
/* Added constraint support that can force a base pair, prohibit a base pair and make single stranded regions */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "data.h"
#include "constants.h"
#include "main-c.h"
#include "algorithms.h"
#ifdef _OPENMP   /* The compiler automatically recognizes openmp support and includes the file accordingly*/
#include "omp.h"
#endif

#define DEBUG 0

#define WM(i,j) WM[j][i]  /* This pragma is defined for readability purpose.*/

unsigned int chPairKey;

//Constraint arrays
int plen = 0, flen = 0, sslen = 0;
int *pbpi, *pbpj, *fbpi, *fbpj, *ss;

/* This function calculates chPairKey to be processed by function chPair. Defined by Professor Bader. */
void init_chPair() {
	int i, j;

	chPairKey = 0;
	for (i = 0; i < 4; i++)
		for (j = 0; j < 4; j++)
			chPairKey += checkPair(i, j) << ((i << 2) + j);
}

int update_chPair(int i, int j) 
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


/* This pragma returns 1 if  base b1 and b2 can pair up, otherwise returns 0, using chPairKey calculated in init_chPair function. Here b1 and b2 are 0-3 to represent one of the four nucleotides A, C, G and U. */
#if 0
inline
int chPair(int b1, int b2) {
	return (chPairKey & (1 << ((b1<<2) + b2)));
}
#else
#define chPair(a, b)  (chPairKey & (1 << (((a)<<2) + (b))))  /* Please try to run this, to understand this statement. Defined by Professor Bader. */
#endif

/* Initialize variables.*/
void initTables(int len) {

	int i, j;
	int LLL;

#if 0
	int z = (len)*(len+1)/2 + 1;

	V = new int[z];
	indx = new int[len+1];
#endif

#if DEBUG
#ifdef DYNALLOC
	if (len != LENGTH-1)
		fprintf(stderr,"ERROR: in initTables, len (%5d) != LENGTH-1 (%5d)\n",len,LENGTH-1);
#endif
#endif

	init_chPair();

	for (i = 0; i < LENGTH; i++) {
		W[i] = INFINITY_; /* Initializing W array with INFINITY make sure that an unfolded sequence will have a large +ve value for free energy - INIFINITY*/
		constraints[i] = 0;
#if 0
		indx[i] = (LENGTH-1)*(i-1) - (i*(i-1))/2;
		indx[i] = (len)*(i-1) - (i*(i-1))/2;
#endif
		for (j = 0; j < LENGTH; j++) {
			VBI[i][j] = INFINITY_;
			VM[i][j] = INFINITY_;
			WM[i][j] = INFINITY_;
		}
	}

	LLL = (LENGTH - 1) * (LENGTH) / 2 + 1;

	for (i = 0; i < LLL; i++)
		V[i] = INFINITY_;

	/*The array V is mapped from 2D to 1D and indexed using the indx array. This mapping helps removing the space wasted for j < i*/
	for (i = 0; i <= LENGTH - 1; i++)
		indx[i] = (len) * (i - 1) - (i * (i - 1)) / 2;

	return;
}

//check if single stranded region is allowed with the given constraints
int checkSS(int i, int j) {

	int it;
	for (it = i + 1; it < j; it++) {
		if (constraints[it] > 0)
			return 1;
	}
	return 0;

}

int calculate(int len, int **forceList, int **prohibitList, int forcelen, int prohibitlen) {
	int b, i, j, it, k;
	
	printf("chPairKey %d ", chPairKey);

	for(i=1;i<=len;i++) 
	{
		if(RNA1[i]=='N') 
			constraints[i] = -1;
	}	

	if (prohibitlen != 0) 
	{
		for (it = 0; it < prohibitlen; it++) 
		{
			for(k= 1; k <= prohibitList[it][2];k++)
			{
				constraints[prohibitList[it][0]+k-1] = -1;
				if(prohibitList[it][1]!=0)
				{
					constraints[prohibitList[it][1]+1-k] = -1;
				}
			}
		}
	}

	if (forcelen != 0) 
	{
		printf("Running with constraints\n");
		for (it = 0; it < forcelen; it++) 
		{
			for(k=1; k <= forceList[it][2];k++)
			{
				if (!chPair(RNA[forceList[it][0]+k-1], RNA[forceList[it][1]-k+1])) 
				{
					printf("Can't constrain (%d,%d)\n", forceList[it][0]+k-1,
							forceList[it][1]-k+1);
					continue;
				}
				constraints[forceList[it][0]+k-1] = forceList[it][1]+1-k;
				constraints[forceList[it][1]+1-k] = forceList[it][0]+k-1;
			}
				//printf("(%d,%d)\n", forceList[it][0], forceList[it][1]);
		}
	}

#if 1
#ifdef _OPENMP
#pragma omp parallel
#pragma omp master
	{
		fprintf(stdout,"\n\n");
		fprintf(stdout,"Running with %3d OpenMP thread",omp_get_num_threads());
		if (omp_get_num_threads()>1) fprintf(stdout,"s");
		fprintf(stdout,".\n\n");
	}
#endif
#endif

    //printf("starting.......\n");

	/* Here b-1 is the length of the segment closed with (i,j) base pair. We assume the minimum size of a hairpin loop closed with (i,j) equal to 3.*/

	/* For b = 4 to 6, hairpin loops and at b = 6 stack loops are possible. So, only WM, and V array are needs to be calculated.
	 * If (i,j) can not pair up then only WM needs to be calculated.
	 * */
	for (b = 4; b <= 6; b++) {
#ifdef _OPENMP
#pragma omp parallel for private (i,j) schedule(guided)
		/* OpenMP syntex to parallelize the for loop. Guided scheduling strategy works better because there may not be equla amount of work for every (i,j)
		 * Please look at the conference paper on GTfold for information regarding how to parallelize the code. Also, note that the for every value of b calculation has to go sequentially. However, the calculation for a perticular value of b, which corresponds to calculating on a single diagonal.
		 * */
#endif
		for (i = 1; i <= len - b; i++) {
			j = i + b;
			//if (constraints[i] == -1 && constraints[j] == -1)
			//	continue;
			if (chPair(RNA[i], RNA[j])) /* Check if bases i and j pair up or not */
				calcVWM(i, j, INFINITY_, INFINITY_); /* Calculates V and WM array for element (i,j)*/
			else
				calcWM(i, j); /* Calculates WM array for element (i,j)*/
		}
	}

	/* Please note that computations of internal loops using speedup algorithm has to be done for every closing base pair (i,j) even if it is not capable of pairing up.
	 * To take care of this, both cases have been separated using a boolean variable ILSA*/

	if (ILSA == FALSE) { /* If we are executing internal loop speedup algorithm (ILSA) */
		/* For b=7 to 10, base pair (i,j) is not able to form multiloops. */
		for (b = 7; b <= 10; b++) {
#ifdef _OPENMP
#pragma omp parallel for private (i,j) schedule(guided)
#endif
			for (i = 1; i <= len - b; i++) {
				j = i + b;
				//if (constraints[i] == -1 && constraints[j] == -1)
				//	continue;
				if (chPair(RNA[i], RNA[j])) {
					calcVBI(i, j); /* Calculates VBI element at (i,j) */
					calcVWM(i, j, VBI[i][j], INFINITY_); /* Calculates V and WM arrays*/
				} else
					calcWM(i, j); /* Calculates WM element at (i,j) */
			}
		}

		for (b = 11; b <= len - 1; b++) {
#ifdef _OPENMP
#pragma omp parallel for private (i,j) schedule(guided)
#endif
			for (i = 1; i <= len - b; i++) {
				j = i + b;
                //printf("%d %d: %d %d\n", i, j, constraints[i], constraints[j]);
				//if (constraints[i] == -1 && constraints[j] == -1)
				//	continue;
				if (chPair(RNA[i], RNA[j])) {
					calcVBIVMVWM(i, j); /* Calculates VBI, VM, V and WM elements at (i,j) */
				} else
					calcWM(i, j); /* Calculates WM element at (i,j) */
			}
		}
	} else { /* If we are executing with ILSA - Internal loop speedup algorithm */
		for (b = 7; b <= 10; b++) {
#ifdef _OPENMP
#pragma omp parallel for private (i,j) schedule(guided)
#endif
			for (i = 1; i <= len - b; i++) {
				j = i + b;
				calcVBIS(i, j); /* Calculates VBI[i][j] array with Internal loop speedup algorithm (ILSA) */
				//if (constraints[i] == -1 && constraints[j] == -1)
				//	continue;
				if (chPair(RNA[i], RNA[j])) {
					calcVWM(i, j, VBI[i][j], INFINITY_); /* Calculates V and WM element at (i,j) */
				} else {
					calcWM(i, j);
				} /* Calculates WM element at (i,j) */
			}
		}

		for (b = 11; b <= len - 1; b++) {
#ifdef _OPENMP
#pragma omp parallel for private (i,j) schedule(guided)
#endif
			for (i = 1; i <= len - b; i++) {
				j = i + b;
				calcVBIS(i, j); /* Calculation of VBI array at (i,j) - Done in both cases whether (i,j) pairs up or not*/
				//if (constraints[i] == -1 && constraints[j] == -1)
				//	continue;
				if (chPair(RNA[i], RNA[j])) {
					calcVMVWM(i, j); /* Calculation of VM, V, WM in Order at (i,j)*/
				} else
					calcWM(i, j); /* Calculation of WM at (i,j)*/
			}
		}
	}

/*
    for(j=2; j<=len; j++){
        for (i=j-1; i>0; i--){
		    printf("%d, (%d,%d), WM: %d\n", j-i, i, j, WM[i][j]);
        }
    }
*/            

	for (j = 5; j <= len; j++) /* Recurssion relation for W array does not depend upon any other array, so can be done after the computation of other arrays are finished.*/
		calcW(j);

	//printV(len);

	return W[len];
}
/* This function calculates the optimal energy of internal loops closed with base pair (i,j) using a heuristic, which limits their size to a constant value - MAXLOOP
 An internal loop contains one closing base pair (i,j) and one enclosed base pair (ip,jp). This function searches for the best enclosed base pair for the closing base pair (i,j) within the given window limited by MAXLOOP
 */
void calcVBI(int i, int j) {

	int ip, jp, temp, VBIij, thres1;

	VBIij = INFINITY_;

	if ((constraints[i] > 0 && constraints[i] != j) || (constraints[j] > 0
			&& constraints[j] != i) || constraints[i] == -1 || constraints[j]
			                                                               == -1)
		return;

	/* Having ip = i+1 and jp = j-1, creates a stack loop. Stack loops are taken care separately in the calculation of V  using eS() function. Therefore, for ip=i+1, the jp value should be lesser than or equal to j-2.*/
	ip = i + 1;
	thres1 = MAX((j - 1) + (ip - i - 1) - MAXLOOP, ip + 4); /* Minimum size of the hairpin loop which the enclosed base pair (ip,jp) can close is 3 that results in the minimum value of jp = ip+4 */
	for (jp = thres1; jp <= j - 2; jp++) {
		//May need to check the constraint condition here
		if (chPair(RNA[ip], RNA[jp])) {
			if (checkSS(i, ip) || checkSS(jp, j))
				continue;
			temp = eL(i, j, ip, jp) + V[indx[ip] + jp]; /* Energy of internal loop closed by (i,j) and (ip,jp) + the optimal energy of the substructure closed by (ip,jp)*/
			if (VBIij > temp)
				VBIij = temp;
		}
	}

	for (ip = i + 2; ip <= i + MAXLOOP + 1; ip++) {
		thres1 = MAX((j - 1) + (ip - i - 1) - MAXLOOP, ip + 4); /* Minimum size of a hairpin loop is 3, so start jp from ip+4*/
		//May need to check for forced constraints here
		for (jp = thres1; jp <= j - 1; jp++) {
			if (chPair(RNA[ip], RNA[jp])) {
				if (checkSS(i, ip) || checkSS(jp, j))
					continue;
				temp = eL(i, j, ip, jp) + V[indx[ip] + jp]; /* Energy of internal loop closed by (i,j) and (ip,jp) + the optimal energy of the substructure closed by (ip,jp)*/
				if (VBIij > temp)
					VBIij = temp;
			}
		}
	}

	VBI[i][j] = VBIij;
	return;
}

/* Amrita: - Internal loop speedup algorithm */
/* Calculation of internal loops using internal loop speedup algorithm. The algorithm calculates the optimal loop closed with base pair (i,j) */
void calcVBIS(int i, int j) {

	int ip, jp, E, VBIij, c = 3, b, len = LENGTH - 1, E1, E2, g; /* ip and jp form enclosed base pairs and c is a small constant currently taken as 3. The loops having one or both sides smaller than c are calculated as special cases. */

	if ((constraints[i] > 0 && constraints[i] != j) || (constraints[j] > 0
			&& constraints[j] != i) || constraints[i] == -1 || constraints[j]
			                                                               == -1)
		return;

	VBIij = VBI[i][j];

	/*Case1: Loops having first side shorter than c and second side has all allowable sizes */
	/* Having ip = i+1 and jp=j-1 creates a stack which is considered separately in eS function. So here the max value of jp could be j-2*/
	ip = i + 1;
	for (jp = ip + 4; jp <= j - 2; jp++) {
		if (chPair(RNA[ip], RNA[jp])) {
			E = eL(i, j, ip, jp) + V[indx[ip] + jp];
			if (VBIij > E)
				VBIij = E;
		}
	}

	for (ip = i + 2; ip <= i + c; ip++) {
		for (jp = ip + 4; jp <= j - 1; jp++) { /* Minimum size of a hairpin loop is 3.*/
			if (chPair(RNA[ip], RNA[jp])) {
				E = eL(i, j, ip, jp) + V[indx[ip] + jp];
				if (VBIij > E)
					VBIij = E;
			}
		}
	}

	/*Case 2: When the first side is greater or equal to c but the second side is smaller */
	for (ip = i + c + 1; ip < j - 1; ip++) {
		for (jp = j - c; jp <= j - 1 && jp >= ip + 4; jp++) { /* Minimum size of a hairpin loop is 3.*/
			if (chPair(RNA[ip], RNA[jp])) {
				E = eL(i, j, ip, jp) + V[indx[ip] + jp];
				if (VBIij > E)
					VBIij = E;
			}
		}
	}

	/* Case 3: General Case - when both sides of internal loops are greater than or equal to c */
	/*Base cases for this (i,j) are g=j-i-2c-3 and j-i-2c-4, gap values should always be greater than or equal to 3*/
	/* First base case - both sides of the loop are  equal to c.*/
	ip = i + c + 1;
	jp = j - c - 1;

	g = jp - ip - 1;
	if (g < 3) {
		VBI[i][j] = VBIij;
		return;
	} /* if g is lesser than 3, then you don't extend this value. In this case the second base case of g-1=j-i-2c-4 will also not make a valid gap value. */

	E = eL(i, j, ip, jp) + V[indx[ip] + jp];
	if (VBIij > E)
		VBIij = E;

	/* Extend this base case for all closing base pairs of the form ( i-b, j+b ) */

	for (b = 1; b <= MIN(i - 1, len - j); b++) {

		E = eL(i - b, j + b, ip, jp) + V[indx[ip] + jp];

		/* Two more options for base pair (i-b,j+b), which are introduced by having one of the side of the resultant internal loop exactly equal to c */
		/* Second side is c - with closing base pair (i-b,j+b)*/
		int ip1 = i + c + 1 + b;
		int jp1 = (j + b) - c - 1;
		E1 = eL(i - b, j + b, ip1, jp1) + V[indx[ip1] + jp1];

		/* First side is c - with closing base pair (i-b,j+b)*/
		int ip2 = (i - b) + c + 1;
		int jp2 = j - c - 1 - b;
		E2 = eL(i - b, j + b, ip2, jp2) + V[indx[ip2] + jp2];

		if (E > E1) {
			E = E1;
			ip = ip1;
			jp = jp1;
		}
		if (E > E2) {
			E = E2;
			ip = ip2;
			jp = jp2;
		}
		if (VBI[i - b][j + b] > E) {
			VBI[i - b][j + b] = E;
		}
	}

	if (g == 3) {
		VBI[i][j] = VBIij;
		return;
	} /* In this case the gap g-1=2, which should not be extended.*/
	ip = i + c + 2;
	jp = j - c - 1;
	E1 = eL(i, j, ip, jp) + V[indx[ip] + jp];
	E2 = eL(i, j, i + c + 1, j - c - 2) + V[indx[i + c + 1] + j - c - 2];
	if (VBIij > E1)
		VBIij = E1;
	if (VBIij > E2)
		VBIij = E2;

	if (E2 < E1) {
		ip = i + c + 1;
		jp = j - c - 2;
	}

	for (b = 1; b <= MIN(i - 1, len - j); b++) {
		E = eL(i - b, j + b, ip, jp) + V[indx[ip] + jp];

		/*Two more options for base pair (i-b,j+b), having one of the sides equal to c */
		/* First side is equal to c*/
		int ip1 = (i - b) + c + 1;
		int jp1 = (j) - c - 2 - b;
		E1 = eL(i - b, j + b, ip1, jp1) + V[indx[ip1] + jp1];

		/* Second side is equal to c*/
		int ip2 = i + c + 2 + b;
		int jp2 = (j + b) - c - 1;
		E2 = eL(i - b, j + b, ip2, jp2) + V[indx[ip2] + jp2];

		if (E > E1) {
			E = E1;
			ip = ip1;
			jp = jp1;
		}
		if (E > E2) {
			E = E2;
			ip = ip2;
			jp = jp2;
		}
		if (VBI[i - b][j + b] > E) {
			VBI[i - b][j + b] = E;
		}
	}

	VBI[i][j] = VBIij;
}

/* Function for calculating the value of WM(i,j)*/
void calcWM(int i, int j) {

	int b = multConst[2], c = multConst[1]; /* b is the branch penalty and c is penalty for single bases for multiloops*/
	int h;
	/* WMidjd = dangling base on both ith and jth side.  WMidj = dangling base on ith side. WMijd = dangling base on jth side. WMij = no dangling base on both sides  */
	int WMidjd, WMidj, WMijd, WMij, WMijp;
	int rnai, rnaj;
	rnai = RNA[i]; /* Read the value of RNA[i] and RNA[j] in register to make the program execute faster.*/
	rnaj = RNA[j];

	WMijp = INFINITY_; /* See, the data flow through the function - how it has been calculated. */

	/* Minimum size of a hairpin loop is 3, that makes the starting limit of h=i+4 and end limit of j-5*/
	for (h = i + 4; h < j - 4; h++) {
		int temp = WM[i][h] + WM(h+1,j);
		if (temp <= WMijp)
			WMijp = temp;
	}

	WMidjd = INFINITY_;
	WMidj = INFINITY_;
	WMijd = INFINITY_;
	WMij = INFINITY_;

	/* If base i and j pair up. */
	WMij = V[indx[i] + j] + auPen(rnai, rnaj) + b;
	/* If base i+1 and j pair up. Add the dangling interaction energy of base pair (i+1,j) with base i being on the 3' end */
	if (constraints[i] <= 0)
		WMidj = V[indx[i + 1] + j] + dangle[rnaj][RNA[i + 1]][rnai][1] + auPen(
				RNA[i + 1], rnaj) + b + c;
	/* If base i and j-1 pair up. Add the dangling interaction energy of base pair (i,j-1) with base j being on the 5' end */
	if (constraints[j] <= 0)
		WMijd = V[indx[i] + j - 1] + dangle[RNA[j - 1]][rnai][rnaj][0] + auPen(
				rnai, RNA[j - 1]) + b + c;
	/* If base i+1 and j-1 pair up. Add the dangling interaction energy of base pair (i+1,j-1) with base i on the 3' and base j on the 5' end.*/
	if (constraints[i] <= 0 && constraints[j] <= 0)
		WMidjd = V[indx[i + 1] + j - 1]
		           + dangle[RNA[j - 1]][RNA[i + 1]][rnai][1]
		                                                  + dangle[RNA[j - 1]][RNA[i + 1]][rnaj][0] + auPen(RNA[i + 1],
		                                                		  RNA[j - 1]) + b + 2* c ;

	//	if(i==6 && j==11)
	//	  printf("(%d,%d), WMij: %d, WMidj: %d, WMijd: %d, WMidjd: %d, constraints: (%d,%d)\n", i, j, WMij, WMidj, WMijd, WMidjd, constraints[i], constraints[j]);

	/* Take the minimum of all of the terms */
	WMij = MIN(MIN(WMij, WMidj), MIN(WMijd, WMidjd));

	int WMsip1j = INFINITY_;
	int WMsijm1 = INFINITY_;

	if (constraints[i] <= 0)
		WMsip1j = WM[i + 1][j];

	if (constraints[j] <= 0)
		WMsijm1 = WM[i][j - 1];


	WMij = MIN(MIN(WMsip1j + c, WMsijm1 + c), WMij);
	WMij = MIN(WMijp, WMij);

	//printf("%d, (%d,%d), WM: %d\n", j-i, i, j, WMij);

	WM[i][j] = WMij;
	WM(i,j) = WMij; /* Extra instruction. NOTE that - by having this instruction we are making WM array symmetric. The macro will convert this instruction into WM[j][i] = WM[i][j] making WM array symmetric.*/
	return;
}

/* Function used for calculating the value of V and WM for a given i,j pair. Calculation of WM at (i,j) requires value of V at (i,j)*/
void calcVWM(int i, int j, int VBIij, int VMij) {
	int a, b, c, h, Vij, eh, es;
	int WMidjd, WMidj, WMijd, WMij, WMijp;
	int rnai, rnaj;

	rnai = RNA[i];
	rnaj = RNA[j];

	WMidjd = INFINITY_;
	WMidj = INFINITY_;
	WMijd = INFINITY_;
	WMij = INFINITY_;

	/* V starts */
	eh = eH(i, j); /* Energy of a hairpin loop */
	es = eS(i, j); /* Energy of a stack, with (i,j) and (i-1,j+1) base pairs.*/
	if (es == 0) {
		es = INFINITY_;
	} else
		es += V[indx[i + 1] + j - 1];

	Vij = MIN(MIN(eh, es), MIN(VBIij, VMij));

	if ((constraints[i] > 0 && constraints[i] != j) || (constraints[j] > 0
			&& constraints[j] != i) || constraints[i] == -1 || constraints[j]
			                                                               == -1)
		Vij = INFINITY_;

	// printf("%d, V(%d,%d): eh: %d, es: %d, vbi: %d, vm: %d, v: %d\n", j-i, i, j, eh, es, VBIij, VMij, Vij);

	V[indx[i] + j] = Vij;

	if (NOISOLATE == TRUE && Vij < INFINITY_) {
		//Check if i+1,j-1 have paired
		if (V[indx[i + 1] + j - 1] > INFINITY_ - SMALLINFTY_) {
			//If not then check for i-1, j+1

			//Isolated base pairs look ahead
			int eHL = eH(i - 1, j + 1);
			int eSL = eS(i - 1, j + 1) + V[indx[i] + j];
			if (i - 1 == 0) {
				eSL = 0;
				eHL = 0;
			}
			int Vijl = (eHL < eSL) ? eHL : eSL;
			// printf("(%d,%d): Hairpin: %d, Stack: %d, Lookahead: %d\n", i, j, eHL, eSL, Vijl);

			if (Vijl > INFINITY_ - SMALLINFTY_)
				//isolated base pair found.. setting energy to infinity
				V[indx[i] + j] = Vij = INFINITY_;
		}
	}

#if DEBUG
	if (indx[i]+j > (LENGTH-1)*(LENGTH)/2)
		fprintf(stderr,"ERROR: in calcVMW: i: %5d  j: %5d\n",i,j);
#endif
	/* V ends */

	/* WM starts */
	a = multConst[0];
	b = multConst[2];
	c = multConst[1];
	WMijp = INFINITY_;

	for (h = i + 4; h <= j - 5; h++) {
		int temp = WM[i][h] + WM(h+1,j);
		if (temp < WMijp)
			WMijp = temp;
	}

	WMij = Vij + auPen(rnai, rnaj) + b;

	if (constraints[i] <= 0)
		WMidj = V[indx[i + 1] + j] + dangle[rnaj][RNA[i + 1]][rnai][1] + auPen(
				RNA[i + 1], rnaj) + b + c;

	if (constraints[j] <= 0)
		WMijd = V[indx[i] + j - 1] + dangle[RNA[j - 1]][rnai][rnaj][0] + auPen(
				rnai, RNA[j - 1]) + b + c;

	if (constraints[i] <= 0 && constraints[j] <= 0)
		WMidjd = V[indx[i + 1] + j - 1]
		           + dangle[RNA[j - 1]][RNA[i + 1]][rnai][1]
		                                                  + dangle[RNA[j - 1]][RNA[i + 1]][rnaj][0] + auPen(RNA[i + 1],
		                                                		  RNA[j - 1]) + b + 2* c ;

	WMij = MIN(MIN(WMij, WMidj), MIN(WMijd, WMidjd));

	int WMsip1j = INFINITY_;
	int WMsijm1 = INFINITY_;

	if (constraints[i] <= 0)
		WMsip1j = WM[i + 1][j];

	if (constraints[j] <= 0)
		WMsijm1 = WM[i][j - 1];

	WMij = MIN(MIN(WMsip1j + c, WMsijm1 + c), WMij);
	WMij = MIN(WMijp, WMij);

	//printf("%d, (%d,%d), WM: %d\n", j-i, i, j, WMij);

	WM[i][j] = WMij;
	WM(i,j) = WMij; /* extra instruction */
	/* WM ends */
	return;
}

/* Function for calculating VM, V, and WM at point (i,j) in the order . Calculation of V at (i,j) requires VBI and VM. Also, calculation of final value of WM(i,j) requires V at (i,j) */

void calcVMVWM(int i, int j) {

	int a = multConst[0] /*offset penalty for multiloops*/,
	b = multConst[2]/*penalty per branch for multiloops*/, c =
		multConst[1]/* Penalty per single base in the multiloop*/,
		a1, h, es;
	int aupen;
	int WMijp, WMidjd, WMidj, WMijd, WMij;
	int VMij, VMijd, VMidj, VMidjd, A_temp;
	int WMip1hm1 /* WM value at i+1 and j-1 */,
	WMip2hm1/* WM value at i+2 and h-1*/,
	WMhjm1/* WM value at h and j-1*/, WMhjm2/* WM value at h and j-2*/,
	WMhp1j /*WM value at h+1 and j*/;
	int rnai, rnaj;
	int tmp1, tmp2;

	rnai = RNA[i];
	rnaj = RNA[j];

	WMidjd = INFINITY_;
	WMidj = INFINITY_;
	WMijd = INFINITY_;
	WMij = INFINITY_;

	/* VM and WM starts */
	aupen = auPen(rnai, rnaj); /* AU or NON GC penalty for base pair (i,j) */
	VMij = INFINITY_;
	VMijd = INFINITY_;
	VMidj = INFINITY_;
	VMidjd = INFINITY_;

	/* Manoj starts */

	WMijp = WM[i][i + 4] + WM(i+5,j);
	a1 = WM[i][i + 5] + WM(i+6,j);
	if (a1 <= WMijp)
		WMijp = a1;

	/* Here we are doing the calculation of VM and WM at point (i,j) concurrently. This for loop calculates the values of VMij, VMidj, VMijd, VMidjd with the value of WMijp. The value of WMijp is needed for calculating the value of WM at (i,j). However, it should be noted that the final value of WM[i][j] requires value of V at (i,j) and will be calculated in the end. */

	/* There are four possibilities for the multiloop closing base pair for the inclusion of dangling energies.
	 * 1) Including the dangling energy of i+1 base and also for base j-1 with the base pair (i,j) closing the multiloop - VMidjd
	 * 2) Including the danlging energy of i+1 base and NOT including the dangling energy of base j-1 with the closing base pair (i,j) - VMidj
	 * 3) NOT including the danlging energy of i+1 base and including the dangling energy of base j-1 with the closing base pair (i,j) - VMijd
	 * 4) NOT including the danlging energy of i+1 base and NOT including the dangling energy of base j-1 with the closing base pair (i,j)-VMij
	 * */
	for (h = i + 6; h <= j - 5; h++) {

		a1 = WM[i][h];
		WMip1hm1 = WM[i + 1][h - 1];
		WMip2hm1 = WM[i + 2][h - 1];
#if 0
		WMhjm1 = WM[h][j-1];
		WMhjm2 = WM[h][j-2];
		WMhp1j = WM[h+1][j];
#else
		WMhjm1 = WM(h,j-1); /* Preprocessor will convert this into WM[j-1][h]. According to the algorithm, the expression should be WM[h][j-1]  -- in which case, as the value of h changes, this will access elements from a column of the WM matrix. To improve run time performance WM array is made symmetric and so WM[h][j-1] = WM[j-1][h] and it will have accesses in a row.*/
		WMhjm2 = WM(h,j-2); /* Same reason as above */
		WMhp1j = WM(h+1,j); /* Same reason as above */
#endif

		/* WM starts */
		a1 += WMhp1j;
		if (a1 <= WMijp)
			WMijp = a1;
		/* WM ends */

		/* Calculation of the four options for VM*/
		A_temp = WMip1hm1 + WMhjm1;
		if ((A_temp <= VMij))
			VMij = A_temp;

		A_temp = WMip2hm1 + WMhjm1;
		if (A_temp <= VMidj && constraints[i + 1] <= 0)
			VMidj = A_temp;

		A_temp = WMip1hm1 + WMhjm2;
		if (A_temp <= VMijd && constraints[j - 1] <= 0)
			VMijd = A_temp;

		A_temp = WMip2hm1 + WMhjm2;
		if (A_temp <= VMidjd && constraints[i + 1] <= 0 && constraints[j - 1]
		                                                               <= 0)
			VMidjd = A_temp;
	}

#if 0
	VMidj += dangle[rnai][rnaj][RNA[i+1]][0];
	VMijd += dangle[rnai][rnaj][RNA[j-1]][1];
	VMidjd += dangle[rnai][rnaj][RNA[i+1]][0] + dangle[rnai][rnaj][RNA[j-1]][1];
#else
	tmp1 = dangle[rnai][rnaj][RNA[i + 1]][0]; /* Dangling energy of base pair (i,j) with single base i+1 at 5' end */
	tmp2 = dangle[rnai][rnaj][RNA[j - 1]][1]; /* Dangling energy of base pair (i,j) with single base j-1 at 3' end */
	VMidj += (tmp1 + c);
	VMidjd += (tmp1 + c);
	VMijd += (tmp2 + c);
	VMidjd += (tmp2 + c);
#endif

	/* Manoj ends */

	VMij = MIN(MIN(VMij, VMidj), MIN(VMijd, VMidjd));
	VMij = VMij + b + a + aupen;

	if ((constraints[i] > 0 && constraints[i] != j) || (constraints[j] > 0
			&& constraints[j] != i) || constraints[i] == -1 || constraints[j]
			                                                               == -1)
		VMij = INFINITY_;

	VM[i][j] = VMij;
	/* VM ends */

	/* V starts */
	es = eS(i, j);
	if (es == 0) {
		es = INFINITY_;
	} else
		es += V[indx[i + 1] + j - 1];

	int Vij;
	Vij = MIN(MIN(eH(i, j), es), MIN(VBI[i][j], VMij));

	if ((constraints[i] > 0 && constraints[i] != j) || (constraints[j] > 0
			&& constraints[j] != i) || constraints[i] == -1 || constraints[j]
			                                                               == -1)
		Vij = INFINITY_;

	V[indx[i] + j] = Vij;

	if (NOISOLATE == TRUE && Vij < INFINITY_) {
		//Check if i+1,j-1 have paired
		if (V[indx[i + 1] + j - 1] > INFINITY_ - SMALLINFTY_) {
			//If not then check for i-1, j+1

			//Isolated base pairs look ahead
			int eHL = eH(i - 1, j + 1);
			int eSL = eS(i - 1, j + 1) + V[indx[i] + j];
			if (i - 1 == 0) {
				eSL = 0;
				eHL = 0;
			}
			int Vijl = (eHL < eSL) ? eHL : eSL;
			//printf("(%d,%d): Hairpin: %d, Stack: %d, Lookahead: %d\n", i, j, eHL, eSL, Vijl);

			if (Vijl > INFINITY_ - SMALLINFTY_)
				//isolated base pair found.. setting energy to infinity
				V[indx[i] + j] = Vij = INFINITY_;
		}
	}

#if DEBUG
	if (indx[i]+j > (LENGTH-1)*(LENGTH)/2)
		fprintf(stderr,"ERROR: in calcVBIVMVWM: i: %5d  j: %5d\n",i,j);
#endif

	/* V ends */

	/* WM starts */

	//Need to take care of these WMs
	WMij = Vij + auPen(rnai, rnaj) + b;

	if (constraints[i] <= 0)
		WMidj = V[indx[i + 1] + j] + dangle[rnaj][RNA[i + 1]][rnai][1] + auPen(
				RNA[i + 1], rnaj) + b + c;

	if (constraints[j] <= 0)
		WMijd = V[indx[i] + j - 1] + dangle[RNA[j - 1]][rnai][rnaj][0] + auPen(
				rnai, RNA[j - 1]) + b + c;

	if (constraints[i] <= 0 && constraints[j] <= 0)
		WMidjd = V[indx[i + 1] + j - 1]
		           + dangle[RNA[j - 1]][RNA[i + 1]][rnai][1]
		                                                  + dangle[RNA[j - 1]][RNA[i + 1]][rnaj][0] + auPen(RNA[i + 1],
		                                                		  RNA[j - 1]) + b + 2* c ;

	WMij = MIN(MIN(WMij, WMidj), MIN(WMijd, WMidjd));

	int WMsip1j = INFINITY_;
	int WMsijm1 = INFINITY_;

	if (constraints[i] <= 0)
		WMsip1j = WM[i + 1][j];

	if (constraints[j] <= 0)
		WMsijm1 = WM[i][j - 1];

	WMij = MIN(MIN(WMsip1j + c, WMsijm1 + c), WMij);

	//printf("%d, (%d,%d), WM: %d\n", j-i, i, j, WMij);

	WM[i][j] = MIN(WMijp, WMij);
	WM(i,j) = WM[i][j]; /* extra instruction */
	/* WM ends */

	return;
}

/* Calculation of VBI, VM, V and WM in the order. Calculation of V at (i,j) requires VBI and VM. Also, calculation of final value of WM(i,j) requires V at (i,j) */
void calcVBIVMVWM(int i, int j) {

	int ip, jp, temp, VBIij, thres1;
	int
	a = multConst[0] /*a is an offset penalty for multiloops */,
	b = multConst[2]/* b is penalty for multiloop branches, one per branch*/,
	c = multConst[1] /*Penalty for single stranded nucleotides in the multiloops*/,
	a1, h, es;
	int aupen;
	int WMijp, WMidjd, WMidj, WMijd, WMij;
	int VMij, VMijd, VMidj, VMidjd, A_temp;
	int WMip1hm1 /* WM value at i+1 and h-1 */,
	WMip2hm1 /* WM value at i+2 and j-1*/,
	WMhjm1 /*WM value at h and j-1*/,
	WMhjm2 /* WM value at h and j-2*/, WMhp1j /*WM value at h+1 and j*/;
	int rnai, rnaj;
	int tmp1, tmp2;

	rnai = RNA[i];
	rnaj = RNA[j];

	WMidjd = INFINITY_;
	WMidj = INFINITY_;
	WMijd = INFINITY_;
	WMij = INFINITY_;

	/* VBI starts */
	/* Look at the calcVBI function for explanation of internal loop calculations*/
	VBIij = INFINITY_;

	int ifinal, jfinal;
	/* ip=i+1, jp=j-1 closes a stack, so we should set the jp limit till j-2, in the following loop*/
	ip = i + 1;
	thres1 = MAX((j - 1) + (ip - i - 1) - MAXLOOP, ip + 4); /* Minimum size of a hairpin loop is 3. So, start jp from ip+4 */
	for (jp = thres1; jp <= j - 2; jp++) {
		if (chPair(RNA[ip], RNA[jp])) {
			if (checkSS(i, ip) || checkSS(jp, j))
				continue;
			temp = eL(i, j, ip, jp) + V[indx[ip] + jp]; /* Energy of internal loop closed by (i,j) and (ip,jp) + the energy of structure closed by (ip, jp)*/
			if (VBIij > temp) {
				VBIij = temp;
				ifinal = ip;
				jfinal = jp;
			}
		}
	}

	for (ip = i + 2; ip <= i + MAXLOOP + 1; ip++) {
		thres1 = MAX((j - 1) + (ip - i - 1) - MAXLOOP, ip + 4);
		for (jp = thres1; jp <= j - 1; jp++) {
			if (chPair(RNA[ip], RNA[jp])) {
				if (checkSS(i, ip) || checkSS(jp, j))
					continue;
				temp = eL(i, j, ip, jp) + V[indx[ip] + jp];
				if (VBIij > temp) {
					VBIij = temp;
					ifinal = ip;
					jfinal = jp;
				}
			}
		}
	}

	if ((constraints[i] > 0 && constraints[i] != j) || (constraints[j] > 0
			&& constraints[j] != i) || constraints[i] == -1 || constraints[j]
			                                                               == -1)
		VBIij = INFINITY_;

	VBI[i][j] = VBIij;
	/* VBI ends */

	/* VM and WM starts */
	aupen = auPen(rnai, rnaj);

	/* Here we are doing the calculation of VM and WM at point (i,j) concurrently. The following for loop calculates the values of VMij, VMidj, VMijd, VMidjd with the value of WMijp. The value of WMijp is needed for calculating the value of WM at (i,j). However, it should be noted that the final value of WM[i][j] requires value of V at (i,j) which itself require VM(i,j), and will be calculated in the end. */

	/* There are four possibilities for the multiloop closing base pair for the inclusion of dangling energies.
	 * 1) Including the dangling energy of i+1 base and also for base j-1 with the base pair (i,j) closing the multiloop - VMidjd
	 * 2) Including the danlging energy of i+1 base and NOT including the dangling energy of base j-1 with the closing base pair (i,j) - VMidj
	 * 3) NOT including the danlging energy of i+1 base and including the dangling energy of base j-1 with the closing base pair (i,j) - VMijd
	 * 4) NOT including the danlging energy of i+1 base and NOT including the dangling energy of base j-1 with the closing base pair (i,j)-VMij
	 * */

	VMij = INFINITY_;
	VMijd = INFINITY_;
	VMidj = INFINITY_;
	VMidjd = INFINITY_;

	/* Manoj starts */
	/* Merged calculations of WM(i,j) and VM(i,j)*/
	WMijp = WM[i][i + 4] + WM(i+5,j);
	a1 = WM[i][i + 5] + WM(i+6,j);
	if (a1 <= WMijp)
		WMijp = a1;

	for (h = i + 6; h <= j - 5; h++) {

		a1 = WM[i][h];
		WMip1hm1 = WM[i + 1][h - 1];
		WMip2hm1 = WM[i + 2][h - 1];
#if 0
		WMhjm1 = WM[h][j-1];
		WMhjm2 = WM[h][j-2];
		WMhp1j = WM[h+1][j];
#else
		WMhjm1 = WM(h,j-1);
		WMhjm2 = WM(h,j-2);
		WMhp1j = WM(h+1,j);
#endif

		/* WM starts */
		a1 += WMhp1j;
		if (a1 <= WMijp)
			WMijp = a1;
		/* WM ends */

		/* Calculation of the four options for VM*/
		A_temp = WMip1hm1 + WMhjm1;
		if ((A_temp <= VMij))
			VMij = A_temp;

		A_temp = WMip2hm1 + WMhjm1;
		if (A_temp <= VMidj && constraints[i + 1] <= 0)
			VMidj = A_temp;

		A_temp = WMip1hm1 + WMhjm2;
		if (A_temp <= VMijd && constraints[j - 1] <= 0)
			VMijd = A_temp;

		A_temp = WMip2hm1 + WMhjm2;
		if (A_temp <= VMidjd && constraints[i + 1] <= 0 && constraints[j - 1]
		                                                               <= 0)
			VMidjd = A_temp;

		// if(i==28 && j==42)
		// 	printf("(%d,%d,%d), VMij: %d, VMidj: %d, VMijd: %d, VMidjd: %d\n", i,h,j, VMij, VMijd, VMidj, VMidjd);

	}

#if 0
	VMidj += dangle[rnai][rnaj][RNA[i+1]][0];
	VMijd += dangle[rnai][rnaj][RNA[j-1]][1];
	VMidjd += dangle[rnai][rnaj][RNA[i+1]][0] + dangle[rnai][rnaj][RNA[j-1]][1];
#else
	tmp1 = dangle[rnai][rnaj][RNA[i + 1]][0];
	tmp2 = dangle[rnai][rnaj][RNA[j - 1]][1];

	VMidj += tmp1;
	VMidjd += tmp1;
	VMijd += tmp2;
	VMidjd += tmp2;
#endif

	//if(i==28 && j==42)
	//	printf("after for: (%d,%d,%d), VMij: %d, VMidj: %d, VMijd: %d, VMidjd: %d\n", i,h,j, VMij, VMijd, VMidj, VMidjd);

	/* Manoj ends */
	VMij = MIN(MIN(VMij, VMidj), MIN(VMijd, VMidjd));
	VMij = VMij + b + a;
	VMij += aupen;

	if ((constraints[i] > 0 && constraints[i] != j) || (constraints[j] > 0
			&& constraints[j] != i) || constraints[i] == -1 || constraints[j]
			                                                               == -1)
		VMij = INFINITY_;

	VM[i][j] = VMij;
	/* VM ends */

	/* V starts */
	es = eS(i, j); /* Energy of stack closed with (i,j) and (i+1,j-1)*/
	if (es == 0) { /* Amrita: I don't know, if this statement is necessary. NOT DONE BY ME. */
		es = INFINITY_;
	} else
		es += V[indx[i + 1] + j - 1];

	int Vij;
	Vij = MIN(MIN(eH(i, j), es), MIN(VBI[i][j], VMij));

	if ((constraints[i] > 0 && constraints[i] != j) || (constraints[j] > 0
			&& constraints[j] != i) || constraints[i] == -1 || constraints[j]
			                                                               == -1)
		Vij = INFINITY_;

	//printf("%d, V(%d,%d): eh: %d, es: %d, vbi: %d, vm: %d, v: %d\n", j-i, i, j, eH(i,j), es, VBI[i][j], VMij, Vij);

	V[indx[i] + j] = Vij;

	if (NOISOLATE == TRUE && Vij < INFINITY_) {
		//Check if i+1,j-1 have paired
		if (V[indx[i + 1] + j - 1] > INFINITY_ - SMALLINFTY_) {
			//If not then check for i-1, j+1

			//Isolated base pairs look ahead
			int eHL = eH(i - 1, j + 1);
			int eSL = eS(i - 1, j + 1) + V[indx[i] + j];
			if (i - 1 == 0) {
				eSL = 0;
				eHL = 0;
			}
			int Vijl = (eHL < eSL) ? eHL : eSL;
			//printf("(%d,%d): Hairpin: %d, Stack: %d, Lookahead: %d\n", i, j, eHL, eSL, Vijl);

			if (Vijl > INFINITY_ - SMALLINFTY_)
				//isolated base pair found.. setting energy to infinity
				V[indx[i] + j] = Vij = INFINITY_;
		}
	}

#if DEBUG
	if (indx[i]+j > (LENGTH-1)*(LENGTH)/2)
		fprintf(stderr,"ERROR: in calcVBIVMVWM: i: %5d  j: %5d\n",i,j);
#endif

	/* V ends */

	/* WM starts */
	/* No dangling base on any of the side */
	WMij = V[indx[i] + j] + aupen + b;
	/* Dangling base i on 3' end of the base pair (i+1,j) */
	if (constraints[i] <= 0)
		WMidj = V[indx[i + 1] + j] + dangle[rnaj][RNA[i + 1]][rnai][1] + auPen(
				RNA[i + 1], rnaj) + b + c;
	/* Dangling base j on 5' end of the base pair (i,j-1)*/
	if (constraints[j] <= 0)
		WMijd = V[indx[i] + j - 1] + dangle[RNA[j - 1]][rnai][rnaj][0] + auPen(
				rnai, RNA[j - 1]) + b + c;
	/* Dangling base i on the 3' end and base j on the 5' end of the base pair (i+1,j-1)*/
	if (constraints[i] <= 0 && constraints[j] <= 0)
		WMidjd = V[indx[i + 1] + j - 1]
		           + dangle[RNA[j - 1]][RNA[i + 1]][rnai][1]
		                                                  + dangle[RNA[j - 1]][RNA[i + 1]][rnaj][0] + auPen(RNA[i + 1],
		                                                		  RNA[j - 1]) + b + 2* c ;

	// if(i==36 && j==50)
	//	  printf("(%d,%d), WMij: %d, WMidj: %d, WMijd: %d, WMidjd: %d\n", i, j, WMij, WMidj, WMijd, WMidjd);

	WMij = MIN(MIN(WMij, WMidj), MIN(WMijd, WMidjd));
	/* Term WM[i+1][j] takes care of the option when base i is neither pairing up nor playing the role in the dangling energy calculation and we have to add penalty 'c' for base i to remain single.
	 * Term WM[i][j-1] takes care of the option when base j is neither pairing up nor playing the role in the dangling energy calculations and add penalty 'c' for base j to remain single.
	 * */

	int WMsip1j = INFINITY_;
	int WMsijm1 = INFINITY_;

	if (constraints[i] <= 0)
		WMsip1j = WM[i + 1][j];

	if (constraints[j] <= 0)
		WMsijm1 = WM[i][j - 1];

	//	  if(i==36 && j==50)
	//		  printf("(%d,%d), WMij: %d, WMsip1j: %d, WMsijm1: %d\n", i, j, WMij, WMsip1j, WMsijm1);


	WMij = MIN(MIN(WMsip1j + c, WMsijm1 + c), WMij);

	//  if(i==35 && j==50)
	//	  printf("(%d,%d), WMij: %d, WMidj: %d, WMijd: %d, WMidjd: %d\n", i, j, WMij, WMidj, WMijd, WMidjd);

	//printf("%d, (%d,%d), WM: %d\n", j-i, i, j, WMij);

	WM[i][j] = MIN(WMijp, WMij);
	WM(i,j) = WM[i][j]; /* extra instruction */
	/* WM ends */

	return;
}

//problems in this function.. i think this should fix it..
/* Function to calculate the value of W[j]. */
void calcW(int j) {

	int i;
	int Wj, Widjd /*Dangling base on both sides*/,
		Wijd/* Dangling base on jth side.*/,
		Widj/* Dangling base on ith side */,
		Wij/* No dangle base on any of the sides */, Wim1 /* Value of W at (i-1). Set to zero if positive*/;
	int rnai, rnaj;
	int must_branch = 0, besti = 0;

	Wj = INFINITY_;

	rnaj = RNA[j];

	for (i = 1; i < j - 3; i++) {

		Wij = Widjd = Wijd = Widj = INFINITY_;

		// printf("i: %d, j: %d\n", i, j);

# if 0
		Wim1=W[i-1];
#endif
#if 1
		Wim1 = MIN(0, W[i - 1]); /* If W[i-1] >=0, this means that there is a branch contained in the sequence from 1 to i-1. Otherwise W[i-1] will be INFINITY. Here Wim1 is defined in this manner, to make the energy of unfolded sequence as INFINITY. */
#endif

		//Wim1 = W[i - 1];

		rnai = RNA[i];

		/* SH: Calculate the energy with no dangle bases. */
		Wij = V[indx[i] + j] + auPen(rnai, rnaj) + Wim1;
		/* Dangle on both sides of the base pair (i+1,j-1). Add the corresponding energy. */
		if (constraints[i] <= 0 && constraints[j] <= 0)
			Widjd = V[indx[i + 1] + j - 1] + auPen(RNA[i + 1], RNA[j - 1])
				+ dangle[RNA[j - 1]][RNA[i + 1]][rnai][1] + dangle[RNA[j
				- 1]][RNA[i + 1]][rnaj][0] + Wim1;
		/* Single base j dangling on the 5' end of base pair (i,j-1) */
		if (constraints[j] <= 0)
			Wijd = V[indx[i] + j - 1] + auPen(rnai, RNA[j - 1]) + dangle[RNA[j
				- 1]][rnai][rnaj][0] + Wim1;
		/* Single base i dangling on the 3' end of base pair (i+1,j)  */
		if (constraints[i] <= 0)
			Widj = V[indx[i + 1] + j] + auPen(RNA[i + 1], rnaj)
				+ dangle[rnaj][RNA[i + 1]][rnai][1] + Wim1;

		int tmpWj = Wj;
		Wj = MIN(MIN(MIN(Wij, Widjd), MIN(Wijd, Widj)), Wj); /* Take the minimum */
		if (tmpWj != Wj) {
			must_branch = 0;
			besti = i;
		}

		// if(i==46 && j==395){
		//      printf("V(%d, %d): %d, WM: %d\n", 46, j, V[indx[46]+j], WM[46][j]);
		//      printf("Wij: %d, Widjd: %d, Wijd: %d, Widj: %d\n\n", Wij, Widjd, Wijd, Widj);
		//}

		if (Wj < INFINITY_) {
			if (Wj == Wij) {
				if (constraints[i] == j) {
					must_branch = 1;
				}
			} else if (Wj == Widjd) {
				if (constraints[i + 1] == j - 1) {
					must_branch = 1;
				}
			} else if (Wj == Wijd) {
				if (constraints[i] == j - 1) {
					must_branch = 1;
				}
			} else {
				if (constraints[i + 1] == j) {
					must_branch = 1;
				}
			}
		}

	}

	//printf("W%dcal: %d, must_branch: %d, best_i: %d", j, Wj, must_branch, besti);

	/* If jth base is not contributing in the energy calculation of W[j] */
	if (!must_branch) {
		if (Wj > W[j - 1])
			Wj = W[j - 1];
	}

	W[j] = Wj;

	//printf(",W%dset: %d\n", j, W[j]);

	//  if(j==11 || j==35 || j==36)
	// printf("\n*****\nMust branch: %d, W%d: %d\n*****\n", must_branch, j, W[j]);
	return;
}

/* For details on calculating energy of different types of Internal loops, please look at the chapter 3 of MS thesis of Mirela Stefania Andronescu.
 * http://www.cs.ubc.ca/grads/resources/thesis/Nov03/Mirela_Andronescu.pdf
 * */
/* Calculates the energy of internal loop with (i,j) as closing base pair and (ip,jp) as enclosed base pair */
int eL(int i, int j, int ip, int jp) {

	int energy;
	int size1, size2, size;
	int loginc; /* SH: Originally unassiged, but needs to be set to 0 so it doesn't throw off later calculations. */
	int lopsided; /* define the asymmetry of an interior loop */

	energy = INFINITY_;
	loginc = 0;

	/*SH: These calculations used to incorrectly be within the bulge loop code, moved out here. */
	size1 = ip - i - 1;
	size2 = j - jp - 1;
	size = size1 + size2;

	if (size1 == 0 || size2 == 0) {
		if (size > 30) {
			/* AM: Does not depend upon i and j and ip and jp - Stacking Energies */
			loginc = (int) floor(prelog * log((double) size / 30.0));
			energy = bulge[30] + eparam[2] + loginc + auPen(RNA[i], RNA[j])
			+ auPen(RNA[ip], RNA[jp]);
		} else if (size <= 30 && size != 1) {
			/* Does not depend upon i and j and ip and jp - Stacking Energies  */
			energy = bulge[size] + eparam[2];
			energy += auPen(RNA[i], RNA[j]) + auPen(RNA[ip], RNA[jp]);
		} else if (size == 1) {
			energy = stack[fourBaseIndex(RNA[i], RNA[j], RNA[ip], RNA[jp])]
			               + bulge[size] + eparam[2]; /* mans */
		}
	} else {
		/* Internal loop */
		lopsided = abs(size1 - size2);

		if (size > 30) {
			loginc = (int) floor(prelog * log((double) size / 30.0));

			/* Please check what should be the difference in the following two options. Is it correct?*/
			if (!((size1 == 1 || size2 == 1) && gail)) { /* normal internal loop with size > 30*/

				energy = tstki[fourBaseIndex(RNA[i], RNA[j], RNA[i + 1], RNA[j
				                                                             - 1])] + tstki[fourBaseIndex(RNA[jp], RNA[ip], RNA[jp
				                                                                                                                + 1], RNA[ip - 1])] + inter[30] + loginc + eparam[3]
				                                                                                                                                                                  + MIN(maxpen, (lopsided * poppen[MIN(2, MIN(size1,
				                                                                                                                                                                		  size2))]));
			} else { /* if size is more than 30 and it is a grossely asymmetric internal loop and gail is not zero*/
				energy
				= tstki[fourBaseIndex(RNA[i], RNA[j], BASE_A, BASE_A)]
				        + tstki[fourBaseIndex(RNA[jp], RNA[ip], BASE_A,
				        		BASE_A)] + inter[30] + loginc
				        		+ eparam[3] + MIN(maxpen, (lopsided
				        				* poppen[MIN(2, MIN(size1, size2))]));
			}
		}
		/* if size is not > 30, we have a looooot of cases... */
		else if (size1 == 2 && size2 == 2) {
			/* 2x2 internal loop */
			energy
			= iloop22[RNA[i]][RNA[ip]][RNA[j]][RNA[jp]][RNA[i + 1]][RNA[i
			                                                            + 2]][RNA[j - 1]][RNA[j - 2]];
		} else if (size1 == 1 && size2 == 2) {
			energy
			= iloop21[RNA[i]][RNA[j]][RNA[i + 1]][RNA[j - 1]][RNA[j - 2]][RNA[ip]][RNA[jp]];
		} else if (size1 == 2 && size2 == 1) {
			/* 1x2 internal loop */
			energy = iloop21[RNA[jp]][RNA[ip]][RNA[j - 1]][RNA[i + 2]][RNA[i
			                                                               + 1]][RNA[j]][RNA[i]];
		} else if (size == 2) {
			/* 1*1 internal loops */
			energy
			= iloop11[RNA[i]][RNA[i + 1]][RNA[ip]][RNA[j]][RNA[j - 1]][RNA[jp]];
		} else if ((size1 == 1 || size2 == 1) && gail) { /* gail = (Grossly Asymmetric Interior Loop Rule) (on/off <-> 1/0)  */
			energy = tstki[fourBaseIndex(RNA[i], RNA[j], BASE_A, BASE_A)]
			               + tstki[fourBaseIndex(RNA[jp], RNA[ip], BASE_A, BASE_A)]
			                       + inter[size] + loginc + eparam[3] + MIN(maxpen, (lopsided
			                    		   * poppen[MIN(2, MIN(size1, size2))]));
		} else { /* General Internal loops */
			energy
			= tstki[fourBaseIndex(RNA[i], RNA[j], RNA[i + 1],
					RNA[j - 1])] + tstki[fourBaseIndex(RNA[jp],
							RNA[ip], RNA[jp + 1], RNA[ip - 1])] + inter[size]
							                                            + loginc + eparam[3] /* AM: I don't understand this eparam value, I think they do not play any role currently. Please look in loader.cc file, for what value have been assinged to various elements of eparam array */
							                                                              + MIN(maxpen,
							                                                            		  (lopsided * poppen[MIN(2, MIN(size1, size2))])); /*  */
		}
	}

	return energy;
}

/* For details on calculating energy of different types of hairpin loops, please look at the chapter 3 of MS thesis of Mirela Stefania Andronescu.
 * http://www.cs.ubc.ca/grads/resources/thesis/Nov03/Mirela_Andronescu.pdf
 * */
/* SH: Function used to calculate energy of a hairpin loop between i & j. */
/* inline */
int eH(int i, int j) {
	/*  Hairpin loop for all the bases between i and j */
	/*  size for size of the loop, energy is the result, loginc is for the extrapolation for loops bigger than 30 */
	int size;
	int loginc;
	int energy = INFINITY_;
	int key, index, count, tlink, kmult;

	size = j - i - 1; /*  size is the number of bases in the loop, when the closing pair is excluded */

	//checking if single stranded region is allowed with given constraints
	if (checkSS(i, j) || (constraints[i] > 0 && constraints[i] != j)
			|| (constraints[j] > 0 && constraints[j] != i))
		return energy;

	/*  look in hairpin, and be careful that there is only 30 values */

	if (size > 30) {
		loginc = (int) ((prelog) * log(((double) size) / 30.0));
		energy = hairpin[30] + loginc + tstkh[fourBaseIndex(RNA[i], RNA[j],
				RNA[i + 1], RNA[j - 1])] + eparam[4]; /* size penalty + terminal mismatch stacking energy*/
	}

	else if (size <= 30 && size > 4) {
		energy = hairpin[size] + tstkh[fourBaseIndex(RNA[i], RNA[j],
				RNA[i + 1], RNA[j - 1])] + eparam[4]; /* size penalty + terminal mismatch stacking energy*/
	}

	else if (size == 4) {
		/*  tetraloop */
		key = 0;
		tlink = 0;
		for (index = 0; index < 6; ++index) {
			switch (RNA[i + index]) {
			case BASE_A:
				kmult = 1;
				break;
			case BASE_C:
				kmult = 2;
				break;
			case BASE_G:
				kmult = 3;
				break;
			case BASE_U:
				kmult = 4;
				break;
			default:
				kmult = 0;
				fprintf(stderr, "ERROR: in tetraloop calculation\n");
			}
			key += kmult * (int) pow(10.0, 5 - index);
		}
		/*  if the sequence is in tloop, we use this value */
		for (count = 1; count < numoftloops && tlink == 0; ++count) {
			if (key == tloop[count][0]) {
				tlink = tloop[count][1];
			}
		}
		energy = tlink + hairpin[size] + tstkh[fourBaseIndex(RNA[i], RNA[j],
				RNA[i + 1], RNA[j - 1])] + eparam[4];
	}

	else if (size == 3) {
		/*  triloop... For the moment, the file triloop.dat is empty */
		/*  else, should have a treatment like the one if size==4 */
		energy = hairpin[size];
		/* AM: Don't include stacking energy terms for triloopls */
		/* + tstkh[RNA[i]][RNA[j]][RNA[i+1]][RNA[j-1]]  */
		/* + eparam[4]; */
		/*  Must be another penalty for terminal AU... Not sure of this */
		energy += auPen(RNA[i], RNA[j]);
	}

	else if (size < 3 && size != 0) {
		/*  no terminal mismatch */
		energy = hairpin[size] + eparam[4];
		if ((RNA[i] == BASE_A && RNA[j] == BASE_U) || (RNA[i] == BASE_U
				&& RNA[j] == BASE_A)) {
			energy += 6; /*  Seems to be a penalty for terminal AU.  *//* Hairpin Loops of size 3 are not allowed, the term hairpin[size] will result in a very large value.  */
		}
	} else if (size == 0)
		return INFINITY_;

	/*  GGG Bonus => GU closure preceded by GG */
	/*  i-2 = i-1 = i = G, and j = U; i < j */
	if (i > 2) {
		if (RNA[i - 2] == BASE_G && RNA[i - 1] == BASE_G && RNA[i] == BASE_G
				&& RNA[j] == BASE_U) {
			energy += gubonus;
			/*  printf ("\n GGG bonus for i %d j %d ", i, j); */
		}
	}

	/*  Poly-C loop => How many C are needed for being a poly-C loop */
	tlink = 1;
	for (index = 1; (index <= size) && (tlink == 1); ++index) {
		if (RNA[i + index] != BASE_C)
			tlink = 0;
	}
	if (tlink == 1) {
		if (size == 3) {
			energy += c3;
		} else {
			energy += cint + size * cslope;
		}
	}

	return energy;
}

/* SH: Function used to calculate energy of stacked pairs (i.j) & (i+1.j-1). */
/* inline */
int eS(int i, int j) {
	int energy;

	if ((constraints[i] > 0 && constraints[i] != j) || (constraints[j] > 0
			&& constraints[j] != i))
		return INFINITY_;

	//if (j - i <= 3)
	//  return INFINITY_;

	/*  not sure about eparam[1], come from MFold.. = 0 */
	energy = stack[fourBaseIndex(RNA[i], RNA[j], RNA[i + 1], RNA[j - 1])]
	               + eparam[1];

	return energy;
}
