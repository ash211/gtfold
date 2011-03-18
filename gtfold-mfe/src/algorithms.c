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
#include <math.h>
#include <time.h>
#include <sys/time.h>
#include <stdlib.h>
#include <assert.h>

#include "constants.h"
#include "utils.h"
#include "energy.h"
#include "global.h"
#include "algorithms.h"

#ifdef _OPENMP   
#include "omp.h"
#endif

double get_seconds() {
	struct timeval tv;
	gettimeofday(&tv, NULL);
	return (double) tv.tv_sec + (double) tv.tv_usec / 1000000.0;
}

int calculate(int len) { 
	int i, j;
	double t1, t2, tint=0, thair=0,tst=0, tml=0, twm=0;

#if 1
#ifdef _OPENMP
#pragma omp parallel
#pragma omp master
		fprintf(stdout,"\n");
		fprintf(stdout,"Running with %3d OpenMP thread",omp_get_num_threads());
		if (omp_get_num_threads()>1) fprintf(stdout,"s");
		fprintf(stdout,".\n");
#endif
#endif
	//for (b = TURN+1; b <= len-1; b++) {
	for (i = len-TURN-1; i >= 1; i--) {
#ifdef _OPENMP
#pragma omp parallel for private (i,j) schedule(guided)
#endif
		//	for (i = 1; i <= len - b; i++) {
		for (j = i+TURN+1; j <= len; j++) {
			//j = i + b;
			if (canPair(RNA[i], RNA[j])) {
				int p, q, h;
				int newV = INFINITY_;
				int VBIij = INFINITY_;

				t1 = get_seconds();
				newV  = MIN(eH(i, j), newV);
				t2 = get_seconds();
				thair += (t2-t1); 

				t1 = get_seconds();
				int stackEnergy = eS(i, j); 
				//newV= MIN(stackEnergy + V(i+1,j-1), newV);
				newV= MIN(stackEnergy + VV1[j-1] , newV);
				t2 = get_seconds();
				tst += (t2-t1); 

				if (j-i > 6) { 
				t1 = get_seconds();
				// Int Loop BEGIN
				for (p = i+1; p <= MIN(j-2-TURN,i+MAXLOOP+1) ; p++) {
					int minq = j-i+p-MAXLOOP-2;
					if (minq < p+1+TURN) minq = p+1+TURN;
					for (q = minq; q < j; q++) {
						if (!canPair(RNA[p], RNA[q])) continue;
						VBIij = MIN(eL(i, j, p, q) + V(p,q), VBIij);
					}
				}

				VBI(i,j) = VBIij;
				// Int Loop END
				t2 = get_seconds();
				tint += (t2-t1); 
				newV = MIN(newV, VBIij);
				}
				if (j-i > 10) {	
					t1 = get_seconds();

					// Multi Loop BEGIN
					int VMij, VMijd, VMidj, VMidjd;
					VMij = VMijd = VMidj = VMidjd = INFINITY_;

					for (h = i+TURN+1; h <= j-1-TURN; h++) { 
						//VMij = MIN(VMij, WM(i+1,h-1) + WM(h,j-1));		
						VMij = MIN(VMij, WMi1[h-1] + WM(h,j-1));		
						//VMidj = MIN(VMidj, WM(i+2,h-1) + WM(h,j-1));	
						VMidj = MIN(VMidj, WMi2[h-1] + WM(h,j-1));	
						//VMijd = MIN(VMijd, WM(i+1,h-1) + WM(h,j-2));	
						VMijd = MIN(VMijd, WMi1[h-1] + WM(h,j-2));	
						//VMidjd = MIN(VMidjd, WM(i+2,h-1) + WM(h,j-2));	
						VMidjd = MIN(VMidjd, WMi2[h-1] + WM(h,j-2));	
					}

					int d3 = Ed3(i,j,j-1);
					int d5 = Ed5(i,j,i+1);

					VMij = MIN(VMij, VMidj + d5 +Ec);
					VMij = MIN(VMij, VMijd + d3 +Ec);
					VMij = MIN(VMij, VMidjd + d5 +  d3+ 2*Ec);
					VMij = VMij + Ea + Eb + auPenalty(i,j);

					VM(i,j) = VMij;
					newV = MIN(newV, VM(i,j));
					// Multi Loop END

					t2 = get_seconds();
					tml += (t2-t1);
				}
				V(i,j) = newV;

				VV[j] = newV;
			}

			if (j-i > 4) {	
				t1 = get_seconds();
				// WM BEGIN
				int h, newWM = INFINITY_;
				//newWM = MIN(V(i,j) + auPenalty(i,j) + Eb, newWM);
				newWM = MIN(VV[j] + auPenalty(i,j) + Eb, newWM);
				//newWM = MIN(V(i+1,j) + Ed3(j,i+1,i) + auPenalty(i+1,j) + Eb + Ec, newWM);
				newWM = MIN(VV1[j] + Ed3(j,i+1,i) + auPenalty(i+1,j) + Eb + Ec, newWM);
				//newWM = MIN(V(i,j-1) + Ed5(j-1,i,j) + auPenalty(i,j-1) + Eb + Ec, newWM);
				newWM = MIN(VV[j-1] + Ed5(j-1,i,j) + auPenalty(i,j-1) + Eb + Ec, newWM);
				//newWM = MIN(V(i+1,j-1) + Ed3(j-1,i+1,i) + Ed5(j-1,i+1,j) + auPenalty(i+1,j-1) + Eb + 2*Ec, newWM) ;
				newWM = MIN(VV1[j-1] + Ed3(j-1,i+1,i) + Ed5(j-1,i+1,j) + auPenalty(i+1,j-1) + Eb + 2*Ec, newWM) ;

				//newWM = MIN(WM(i+1,j) + Ec, newWM);
				newWM = MIN(WMi1[j] + Ec, newWM);
				//newWM = MIN(WM(i,j-1) + Ec, newWM);
				newWM = MIN(WMi[j-1] + Ec, newWM);

				for (h = i ; h <= j-TURN-1; h++) {
					//newWM = MIN(newWM, WM(i,h) + WM(h+1,j));
					newWM = MIN(newWM, WMi[h] + WM(h+1,j));
				}
				WM(i,j) = WMi[j]=newWM;
				// WM END

				t2 = get_seconds();
				twm += (t2-t1);
			}
		}
		int *FF;
		FF = VV1; VV1=VV; VV=FF;
		FF = WMi2; WMi2 = WMi1; WMi1 = WMi; WMi = FF;
		for (j=1; j<=len; j++) {VV[j]=INFINITY_;WMi[j]=INFINITY_; }
	}

	printf("Total IntLoop time = %fs \n",tint);
	printf("Total Stack time = %fs \n",tst);
	printf("Total HairPin time = %fs \n",thair);
	printf("Total MultiLoop time = %fs \n",tml);
	printf("Total WM time = %fs \n",twm);

	for (j = TURN+2; j <= len; j++)	{
		int i, Wj, Widjd, Wijd, Widj, Wij, Wim1;
		Wj = INFINITY_;

		for (i = 1; i < j-TURN; i++) {
			Wij = Widjd = Wijd = Widj = INFINITY_;
			Wim1 = MIN(0, W[i-1]); 
			Wij = V(i, j) + auPenalty(i, j) + Wim1;
			Widjd = V(i+1,j-1) + auPenalty(i+1,j-1) + Ed3(j-1,i + 1,i) + Ed5(j-1,i+1,j) + Wim1;
			Wijd = V(i,j-1) + auPenalty(i,j- 1) + Ed5(j-1,i,j) + Wim1;
			Widj = V(i + 1, j) + auPenalty(i+1,j) + Ed3(j,i + 1,i) + Wim1;
			Wj = MIN(MIN(MIN(Wij, Widjd), MIN(Wijd, Widj)), Wj); 
		}

		W[j] = MIN(Wj, W[j-1]);
	}

	return W[len];
}
