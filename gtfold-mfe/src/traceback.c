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

/* AUTHOR: Amrita Mathuriya August 2007-January 2009.
 * This file tracebacks only one optimal secondary structure. Note that there may be more than one optimal secondary structures and the generated structure depends upon the choices taken during traceback.
 * This file contains recursive functions. structure is a global integer array to store the optimal structure generated during traceback.
 * Please note that - To understand the functions defined in this file, you first need to understand functions defined in algorithms.c file. As all those functions have a corresponding function here, for traceback the calculations done there.
 * */

/* Modified by Sainath Mallidi August 2009 -  "*/
/* Modified traceback to handle constraints */

/* Modified by Prashant Gaurav <pgaurav3@gatech.edu>, 09/07/2010 */
/* Fixed the incorrect reporting of multiloop energy */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "data.h"
#include "constants.h"
#include "algorithms.h"
#include "traceback.h"

int total_en = 0;

void trace(int len) {
	int i;

	for (i = 0; i < LENGTH; i++)
		structure[i] = 0;

	if (W[len] >= MAXENG) {
		printf("\n No Structure ");
		return;
	}

	traceW(len);
	
	printf("SUM of energy of loops = %12.2f\n", total_en/100.0);
	return;
}

/* Traces W[j] */
void traceW(int j) {

	//printf("Call made with j=%d\n", j);

	int done, i, Wj,Wj_temp;
	int wim1, flag, Widjd, Wijd, Widj, Wij;

	flag = 1;

	done = 0; /* the done variable makes it sure that we are tracebacking the first optimal possibility */
	Wj = INFINITY_;
	int min_i=1;
	if (j == 0 || j == 1) {
		/* W[j] = 0; */
	} else {
		for (i = 1; i < j && !done; i++) {
			//wim1 = W[i - 1];
		    wim1 = MIN ( 0, W[i-1] );

			
			flag = 1;
			if (wim1 != W[i-1])
			    flag = 0;

            //printf("i: %d, j: %d, flag: %d\n", i, j, flag);
			 

			Widjd = INFINITY_;
			Wijd = INFINITY_;
			Widj = INFINITY_;

			Wij = V[indx[i] + j] + auPen(RNA[i], RNA[j]) + wim1;

			if (constraints[i] <= 0 && constraints[j] <= 0)
				Widjd = V[indx[i + 1] + j - 1] + auPen(RNA[i + 1], RNA[j - 1])
				+ dangle[RNA[j - 1]][RNA[i + 1]][RNA[i]][1]
				                                         + dangle[RNA[j - 1]][RNA[i + 1]][RNA[j]][0] + wim1;
			if (constraints[j] <= 0)
				Wijd = V[indx[i] + j - 1] + auPen(RNA[i], RNA[j - 1])
				+ dangle[RNA[j - 1]][RNA[i]][RNA[j]][0] + wim1;
			if (constraints[i] <= 0)
				Widj = V[indx[i + 1] + j] + auPen(RNA[i + 1], RNA[j])
				+ dangle[RNA[j]][RNA[i + 1]][RNA[i]][1] + wim1;
			Wj_temp=Wj;

			Wj = MIN(MIN(MIN(Wij, Widjd), MIN(Wijd, Widj)), Wj);
			if (Wj_temp>Wj) min_i=i;	

			//if(W[j]==Wj)
			// printf("Energies at (%d,%d), W[j]: %d, Wj: %d, Wij: %d, Wijd: %d, Widj: %d, Widjd: %d\n", i, j, W[j], Wj, Wij, Wijd, Widj, Widjd);

			if (W[j] == Wj) {
				if (W[j] == Wij) { /* If the optimal secondary structure contain base pair (i,j) as paired.*/
					done = 1;
					structure[i] = j;
					structure[j] = i;
				    //printf("before V1, i: %d, j: %d\n", i, j);
					traceV(i, j);
					printf("AU Penalty: %12.2f\n",auPen(RNA[i], RNA[j])/100.00);

					if (flag || checkSS(1,i)) // Added this condition because if constraints are present, flag should be overridden
						traceW(i - 1);
					break;
				} else if (W[j] == Widjd && constraints[i] <= 0
						&& constraints[j] <= 0) { /* If base pair (i+1,j-1) is pairing and there is a dangling base on both its sides */
					done = 1;
					structure[i + 1] = j - 1;
					structure[j - 1] = i + 1;
					//printf("before V2, i: %d, j: %d, flag: %d\n", i, j, flag);
					traceV(i + 1, j - 1);
					printf("AU Penalty: %12.2f\nEnergy of dangling base on both the sides: %12.2f  %12.2f\n",auPen(RNA[i + 1], RNA[j - 1])/100.00, dangle[RNA[j - 1]][RNA[i + 1]][RNA[i]][1]/100.00, dangle[RNA[j - 1]][RNA[i + 1]][RNA[j]][0]/100.00);

					if (flag || checkSS(1,i))
						traceW(i - 1);
					break;
				} else if (W[j] == Wijd && constraints[j] <= 0) { /* If base pair (i,j-1) pairs and base j is single stranded. */
					done = 1;
					structure[i] = j - 1;
					structure[j - 1] = i;
					//printf("before V3, i: %d, j: %d\n", i, j);
					traceV(i, j - 1);
					printf("AU Penalty: %12.2f\nEnergy of dangling base: %12.2f\n",auPen(RNA[i], RNA[j - 1])/100.00,dangle[RNA[j - 1]][RNA[i]][RNA[j]][0]/100.00);

					if (flag || checkSS(1,i) )
						traceW(i - 1);
					break;
				} else if (W[j] == Widj && constraints[i] <= 0) { /* If base pair (i+1,j) pairs and base i is single stranded. */
					done = 1;
					structure[i + 1] = j;
					structure[j] = i + 1;
					//printf("before V4, i: %d, j: %d\n", i, j);
					traceV(i + 1, j);
					printf("AU Penalty: %12.2f\nEnergy of dangling base:  %12.2f\n",auPen(RNA[i + 1], RNA[j])/100.00, dangle[RNA[j]][RNA[i + 1]][RNA[i]][1]/100.00);
					if (flag || checkSS(1,i))
						traceW(i - 1);
					break;
				}
			}
		}
		
		if (W[j] == W[j - 1] && !done) /* If jth base is single stranded and also does not have dangling energy contribution */
			traceW(j - 1);

	}

	return;
}


/* Trace the structure inside V[i][j]. This function traces "which type of loop (i,j) base pair is closing" */
int traceV(int i, int j) {

	int a, b, c, d, Vij;

	a = eH(i, j);
	b = eS(i, j) + V[indx[i + 1] + j - 1];
	if (eS(i, j) == 0)
		b = INFINITY_;
	c = VBI[i][j];
	d = VM[i][j];
	Vij = MIN(MIN(a, b), MIN(c, d));

	if (Vij == a && Vij != b && Vij != c && Vij != d) { /* If () a hairpin loop */
		printf("i %5d j %5d Hairpin  Loop %12.2f\n", i, j, eH(i, j)/100.00);
		total_en += eH(i,j);
		return Vij;
	} else if (Vij == b /*&& Vij != a && Vij != c && Vij != d*/) { /* If it forms a stack */
		printf("i %5d j %5d Stack 	 Loop %12.2f\n", i, j, eS(i, j)/100.00);
		total_en += eS(i,j);
		structure[i + 1] = j - 1;
		structure[j - 1] = i + 1;
		traceV(i + 1, j - 1);
		return Vij;
	} else if (Vij == c /*&& Vij != a && Vij != b && Vij != d*/) { /* If it forms an internal loop */
		printf("i %5d j %5d Internal Loop", i, j);
		traceVBI(i, j);
		return Vij;
	} else if (Vij == d && Vij != a && Vij != b && Vij != c) { /* If it forms a multiloop */
		int eVM = traceVM(i, j);
		printf("i %5d j %5d Multi    Loop %12.2f\n", i, j, (Vij-eVM)/100.0);
		total_en += (Vij-eVM);
		return Vij;
	} else { /* It should not be executed */
		if (Vij == c && Vij == d)
			printf("Internal Loop and Multiloops are possible\n");
		printf("i %5d j %5d  More than one structure is possible\n", i, j);
		return 0;
	}
}

/* Traces VBI[i][j] */
int traceVBI(int i, int j) {
	
	int VBIij_temp;
	int ip, jp, el, v;
	int ifinal, jfinal;

	ifinal = 0;
	jfinal = 0;

	for (ip = i + 1; ip < j - 1; ip++) {
		for (jp = ip + 1; jp < j; jp++) { /* Search which internal loop (ip,jp) is closing */
			if (checkSS(i, ip) || checkSS(jp, j))
				continue;
			el = eL(i, j, ip, jp);
			v = V[indx[ip] + jp];
			VBIij_temp = el + v;
			if (VBIij_temp == VBI[i][j]) {
				ifinal = ip;
				jfinal = jp;
				break;
			}
		}
		if (jp != j)
			break;
	}

	structure[ifinal] = jfinal;
	structure[jfinal] = ifinal;
	printf(" %12.2f\n", eL(i, j, ifinal, jfinal)/100.00);
	total_en += eL(i, j, ifinal, jfinal);

	int eVI = traceV(ifinal, jfinal);
	return eVI ;
}

/* Tracing VM[i][j] */
int traceVM(int i, int j) {

	int done;
	int h;
	int a, b;
	int A_temp;
	int eVM = 0;

	done = 0;
	a = eparam[5];
	b = eparam[10]; /* efn2b */

	//printf("  %12d\n", VM[i][j]);

	/*  WM (i+1, j-1) */
	for (h = i + 2; h <= j - 1 && !done; h++) {
		A_temp = WM[i + 1][h - 1] + WM[h][j - 1] + a + b + auPen(RNA[i],
				RNA[j]);
		if (A_temp == VM[i][j]) { /* No dangling bases on any of the sides of base pair (i,j) in the multiloop */
			done = 1;
			eVM += traceWM(i + 1, h - 1);
			eVM += traceWM(h, j - 1);
			break;
		}
	}

	/* WM[i+2, j-1] */
	if (constraints[i + 1] <= 0) {
		for (h = i + 3; h <= j - 1 && !done; h++) {
			A_temp = WM[i + 2][h - 1] + WM[h][j - 1] + a + b + auPen(RNA[i],
					RNA[j]) + dangle[RNA[i]][RNA[j]][RNA[i + 1]][0]; /* IF base i+1 is dangling on the 5' end of the base pair (i,j) */
			if (A_temp == VM[i][j]) {
				done = 1;
				eVM += traceWM(i + 2, h - 1);
				eVM += traceWM(h, j - 1);
				break;
			}
		}
	}

	/* WM[i+1][j-2] */
	if (constraints[j - 1] <= 0) {
		for (h = i + 2; h <= j - 2 && !done; h++) { /* If base j-1 is dangling on 3' end of the base pair (i,j) */
			A_temp = WM[i + 1][h - 1] + WM[h][j - 2] + a + b + auPen(RNA[i],
					RNA[j]) + dangle[RNA[i]][RNA[j]][RNA[j - 1]][1];
			if (A_temp == VM[i][j]) {
				done = 1;
				eVM += traceWM(i + 1, h - 1);
				eVM += traceWM(h, j - 2);
				break;
			}
		}
	}

	/* WM[i+2][j-2] */
	if (constraints[i + 1] <= 0 && constraints[j - 1] <= 0) {
		for (h = i + 3; h <= j - 2 && !done; h++) { /* If base pair (i,j) has dangling bases on both sides. */
			A_temp = WM[i + 2][h - 1] + WM[h][j - 2] + a + b + auPen(RNA[i],
					RNA[j]) + dangle[RNA[i]][RNA[j]][RNA[i + 1]][0]
					                                             + dangle[RNA[i]][RNA[j]][RNA[j - 1]][1];
			if (A_temp == VM[i][j]) {
				done = 1;
				eVM += traceWM(i + 2, h - 1);
				eVM += traceWM(h, j - 2);
				break;
			}
		}
	}

	return eVM;
}

/* Tracing WM[i][j] */
int traceWM(int i, int j) {

	int done;
	int a, b, c;
	int h1, h;
	int eWM = 0; 

	done = 0;

	b = eparam[10]; /* efn2b */
	c = eparam[6]; /* efn2c */

	h1 = 0;

	if (i >= j)
		return 0;
	else {
		for (h = i; h < j && !done; h++) {
			a = WM[i][h] + WM[h + 1][j]; /* If WM(i,j) came from the summation of two WM terms */
			if (a == WM[i][j]) {
				done = 1;
				h1 = h;
				break;
			}
		}
		if (h1 != 0) {
			eWM += traceWM(i, h);
			eWM += traceWM(h + 1, j);
		} else {
			if (WM[i][j] == V[indx[i] + j] + auPen(RNA[i], RNA[j]) + b) { /* If base pair (i,j) pairs*/
				done = 1;
				structure[i] = j;
				structure[j] = i;
				//printf("%12.2f %12.2f %12.2f \n",V[indx[i] + j] , auPen(RNA[i], RNA[j]) , b);
				eWM += traceV(i, j);
			} else if (WM[i][j] == V[indx[i + 1] + j] + dangle[RNA[j]][RNA[i
			                                                               + 1]][RNA[i]][1] + auPen(RNA[i + 1], RNA[j]) + b + c
			                                                               && constraints[i] <= 0) { /* If base pair (i+1, j) forms and base i remain as a dangling base*/
				done = 1;
				//printf("%12.2f %12.2f %12.2f %12.2f %12.2f\n",V[indx[i + 1] + j],dangle[RNA[j]][RNA[i+ 1]][RNA[i]][1] , auPen(RNA[i + 1], RNA[j]) , b , c);
				eWM += traceV(i + 1, j);
				structure[i + 1] = j;
				structure[j] = i + 1;
			} else if (WM[i][j] == V[indx[i] + j - 1]
			                         + dangle[RNA[j - 1]][RNA[i]][RNA[j]][0] + auPen(RNA[i],
			                        		 RNA[j - 1]) + b + c && constraints[j] <= 0) { /* If base pair (i,j-1) forms and base j remain as a dangling base */
				done = 1;
				eWM += traceV(i, j - 1);
				structure[i] = j - 1;
				structure[j - 1] = i;
			} else if (WM[i][j] == V[indx[i + 1] + j - 1]
			                         + dangle[RNA[j - 1]][RNA[i + 1]][RNA[i]][1] + dangle[RNA[j
			                                                                                  - 1]][RNA[i + 1]][RNA[j]][0]
			                                                                                                            + auPen(RNA[i + 1], RNA[j - 1]) + b + 2* c
			                                                                                                            && constraints[i] <= 0 && constraints[j] <= 0) { /* IF base pair (i+1,j-1) forms and base i and base j remain as dangling bases on each side.*/
				done = 1;
//printf("%12.2f %12.2f %12.2f %12.2f %12.2f %12.2f \n",V[indx[i + 1] + j - 1], dangle[RNA[j - 1]][RNA[i + 1]][RNA[i]][1],dangle[RNA[j- 1]][RNA[i + 1]][RNA[j]][0], auPen(RNA[i + 1], RNA[j - 1]) , b , 2* c);
				eWM += traceV(i + 1, j - 1);
				structure[i + 1] = j - 1;
				structure[j - 1] = i + 1;
			} else if (WM[i][j] == WM[i + 1][j] + c && constraints[i] <= 0) { /* If base i does not make energy contribution */
				done = 1;
				eWM += traceWM(i + 1, j);
			} else if (WM[i][j] == WM[i][j - 1] + c && constraints[j] <= 0) { /* If base j does not make energy contribution */
				done = 1;
				eWM += traceWM(i, j - 1);
			}
		}
	}
	return eWM;
}

