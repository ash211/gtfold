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

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include "data.h"
#include "constants.h"
#include "energy.h"
#include "global.h"
#include "traceback.h"
#include "utils.h"


int total_en = 0;

void trace(int len) {
	int i;
	
	for (i = 0; i < len+1; i++)
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
void traceW(int j) 
{
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
			wim1 = MIN(0, W[i-1]);
			flag = 1;
			if (wim1 != W[i-1]) flag = 0;

			Widjd = INFINITY_;
			Wijd = INFINITY_;
			Widj = INFINITY_;
			Wij = V(i,j) + auPenalty(i, j) + wim1;
			Widjd = V(i+1,j-1) + auPenalty(i+1, j-1) + Ed3(j-1,i+1,i) + Ed5(j-1,i+1,j) + wim1;
			Wijd = V(i,j-1) + auPenalty(i,j-1) + Ed5(j-1,i,j) + wim1;
			Widj = V(i+1,j) + auPenalty(i+1,j) + Ed3(j,i+1,i) + wim1;
			Wj_temp=Wj;

			Wj = MIN(MIN(MIN(Wij, Widjd), MIN(Wijd, Widj)), Wj);
			if (Wj_temp>Wj) min_i=i;	

			if (W[j] == Wj) {
				if (W[j] == Wij) { /* If the optimal secondary structure contain base pair (i,j) as paired.*/
					done = 1;
					structure[i] = j;
					structure[j] = i;
					//printf("AU Penalty: %12.2f\n",auPen(RNA[i], RNA[j])/100.00);
					traceV(i, j);

					if (flag ) 
						traceW(i - 1);
					break;
				} else if (W[j] == Widjd) { /* If base pair (i+1,j-1) is pairing and there is a dangling base on both its sides */
					done = 1;
					structure[i + 1] = j - 1;
					structure[j - 1] = i + 1;
					//printf("AU Penalty: %12.2f\nEnergy of dangling base on both the sides: %12.2f  %12.2f\n",
					//auPen(RNA[i + 1], RNA[j - 1])/100.00, dangle[RNA[j - 1]][RNA[i + 1]][RNA[i]][1]/100.00, dangle[RNA[j - 1]][RNA[i + 1]][RNA[j]][0]/100.00);
					traceV(i + 1, j - 1);

					if (flag)// || checkSS(1,i))
						traceW(i - 1);
					break;
				} else if (W[j] == Wijd) { /* If base pair (i,j-1) pairs and base j is single stranded. */
					done = 1;
					structure[i] = j - 1;
					structure[j - 1] = i;
					//	printf("AU Penalty: %12.2f \nEnergy of dangling base: %12.2f \n",
					//	auPen(RNA[i], RNA[j - 1])/100.00,dangle[RNA[j - 1]][RNA[i]][RNA[j]][0]/100.00);
					traceV(i, j - 1);

					if (flag) // || checkSS(1,i) )
						traceW(i - 1);
					break;
				} else if (W[j] == Widj) { /* If base pair (i+1,j) pairs and base i is single stranded. */
					done = 1;
					structure[i + 1] = j;
					structure[j] = i + 1;
					//	printf("AU Penalty: %12.2f\nEnergy of dangling base:  %12.2f\n",
					//	auPen(RNA[i + 1], RNA[j])/100.00, dangle[RNA[j]][RNA[i + 1]][RNA[i]][1]/100.00);
					traceV(i + 1, j);
					if (flag) // || checkSS(1,i))
						traceW(i - 1);
					break;
				}
			}
		}
		
		if (W[j] == W[j - 1] && !done) 
			traceW(j-1);

	}

	return;
}


/* Trace the structure inside V[i][j]. This function traces "which type of loop (i,j) base pair is closing" */
int traceV(int i, int j) {

	int a, b, c, d, Vij;

	a = eH(i, j);
	b = eS(i, j) + V[indx[i + 1] + j - 1];
	if (eS(i, j) == 0) b = INFINITY_;
	c = VBI[i][j];
	d = VM[i][j];
	
	Vij = MIN(MIN(a, b), MIN(c, d));
	
	if (Vij == a && Vij != b && Vij != c && Vij != d) { /* If () a hairpin loop */
		//printf("i %5d j %5d Hairpin  Loop %12.2f\n", i, j, eH(i, j)/100.00);
		total_en += eH(i,j);
		return Vij;
	} else if (Vij == b /*&& Vij != a && Vij != c && Vij != d*/) { /* If it forms a stack */
		//printf("i %5d j %5d Stack 	 Loop %12.2f\n", i, j, eS(i, j)/100.00);
		total_en += eS(i,j);
		structure[i + 1] = j - 1;
		structure[j - 1] = i + 1;
		traceV(i + 1, j - 1);
		return Vij;
	} else if (Vij == c /*&& Vij != a && Vij != b && Vij != d*/) { /* If it forms an internal loop */
		//printf("i %5d j %5d Internal Loop", i, j);
		traceVBI(i, j);
		return Vij;
	} else if (Vij == d && Vij != a && Vij != b && Vij != c) { /* If it forms a multiloop */
		int eVM = traceVM(i, j);
		//printf("i %5d j %5d Multi    Loop %12.2f\n", i, j, (Vij-eVM)/100.0);
		total_en += (Vij-eVM);
		return Vij;
	} 
	return 0;
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
	//printf(" %12.2f\n", eL(i, j, ifinal, jfinal)/100.00);
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

	int VMij = VM[i][j];

	for (h = i + 2; h <= j - 1 && !done; h++) {
		A_temp = WM[i + 1][h - 1] + WM[h][j - 1] + a + b + auPen(RNA[i],
				RNA[j]);
		if (A_temp == VMij) { /* No dangling bases on any of the sides of base pair (i,j) in the multiloop */
			done = 1;
			eVM += traceWM(i + 1, h - 1);
			eVM += traceWM(h, j - 1);
			break;
		}
	}

	for (h = i + 3; h <= j - 1 && !done; h++) {
		A_temp = WM[i + 2][h - 1] + WM[h][j - 1] + a + b + auPenalty(i,j) + Ed5(i,j,i + 1); 
		if (A_temp == VMij) {
			done = 1;
			eVM += traceWM(i + 2, h - 1);
			eVM += traceWM(h, j - 1);
			break;
		}
	}

	for (h = i + 2; h <= j - 2 && !done; h++) { /* If base j-1 is dangling on 3' end of the base pair (i,j) */
		A_temp = WM[i + 1][h - 1] + WM[h][j - 2] + a + b + auPenalty(i, j) + Ed3(i,j,j - 1);
		if (A_temp == VMij) {
			done = 1;
			eVM += traceWM(i + 1, h - 1);
			eVM += traceWM(h, j - 2);
			break;
		}
	}

	for (h = i + 3; h <= j - 2 && !done; h++) { /* If base pair (i,j) has dangling bases on both sides. */
		A_temp = WM[i + 2][h - 1] + WM[h][j - 2] + a + b + auPenalty(i,j) + Ed5(i,j,i + 1) + Ed3(i,j,j - 1);
		if (A_temp == VMij) {
			done = 1;
			eVM += traceWM(i + 2, h - 1);
			eVM += traceWM(h, j - 2);
			break;
		}
	}

	return eVM;
}

/* Tracing WM[i][j] */
int traceWM(int i, int j) {

	int done;
	int h1, h;
	int eWM = 0; 

	done = 0;
	h1 = 0;

	if (i >= j)
		return 0;
	else {
		for (h = i; h < j && !done; h++) {
			int aa = WM[i][h] + WM[h + 1][j]; /* If WM(i,j) came from the summation of two WM terms */
			if (aa == WM[i][j]) {
				done = 1;
				h1 = h;
				break;
			}
		}
		if (h1 != 0) {
			eWM += traceWM(i, h);
			eWM += traceWM(h + 1, j);
		} else {
			if (WM[i][j] == V(i,j) + auPenalty(i, j) + Eb()) { /* If base pair (i,j) pairs*/
				done = 1;
				structure[i] = j;
				structure[j] = i;
				eWM += traceV(i, j);
			} else if (WM[i][j] == V(i+1, j) + Ed3(j,i + 1,i) + auPenalty(i+1, j) + Eb() + Ec()) { 
				done = 1;
				eWM += traceV(i + 1, j);
				structure[i + 1] = j;
				structure[j] = i + 1;
			} else if (WM[i][j] == V(i,j-1) + Ed5(j-1,i,j) + auPenalty(i,j-1) +  Eb() + Ec() ) { 
				done = 1;
				eWM += traceV(i, j - 1);
				structure[i] = j - 1;
				structure[j - 1] = i;
			} else if (WM[i][j] == V(i+1,j-1) + Ed3(j-1,i+1,i) + Ed5(j-1,i+1,j) + auPenalty(i+1, j-1) + Eb() + 2*Ec()) { 
				done = 1;
				eWM += traceV(i + 1, j - 1);
				structure[i + 1] = j - 1;
				structure[j - 1] = i + 1;
			} else if (WM[i][j] == WM[i + 1][j] + Ec() ) { 
				done = 1;
				eWM += traceWM(i + 1, j);
			} else if (WM[i][j] == WM[i][j - 1] + Ec() ) { 
				done = 1;
				eWM += traceWM(i, j - 1);
			}
		}
	}
	return eWM;
}
