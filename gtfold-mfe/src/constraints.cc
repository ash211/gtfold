#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <iostream>
#include <fstream>

#include "global.h"
#include "options.h"
#include "constraints.h"

int* BP;
int** PBP;
int** FBP;

int nPBP;
int nFBP;

static int load_constraints(const char* constr_file, int verbose=0) {
	std::ifstream cfcons;
    fprintf(stdout, "- Running with constraints\n");

    cfcons.open(constr_file, std::ios::in);
    if (cfcons == 0) {
        fprintf(stderr, "Error opening constraint file\n\n");
        cfcons.close();
        return -1;
    }

    char cons[100];

    while (!cfcons.eof()) {
        cfcons.getline(cons, 100);
        if (cons[0] == 'F' || cons[0] == 'f') nFBP++;
        if (cons[0] == 'P' || cons[0] == 'p') nPBP++;
    }
    cfcons.close();

    //fprintf(stdout, "Number of Constraints given: %d\n\n", nFBP + nPBP);
    if (nFBP + nPBP == 0) {
        fprintf(stderr, "No Constraints found.\n\n");
        return -1;
    }

    FBP = (int**) malloc(nFBP*sizeof(int*));
    PBP = (int**) malloc(nPBP*sizeof(int*));

    int fit = 0, pit = 0, it = 0;

    for (it = 0; it < nFBP; it++) {
        FBP[it] = (int*) malloc(2* sizeof (int));
    }
    for(it=0; it < nPBP; it++) {
        PBP[it] = (int*)malloc(2*sizeof(int));
    }
    cfcons.open(constr_file, std::ios::in);

    while(!cfcons.eof()) {
        cfcons.getline(cons,100);
        char *p=strtok(cons, " ");
        p = strtok(0, " ");
        if(cons[0]=='F' || cons[0]=='f') {
            int fit1=0;
            while(p!=0) {
                FBP[fit][fit1++] = atoi(p);
                p = strtok(0, " ");
            }
            fit++;
        }
        if(cons[0]=='P' || cons[0]=='p') {
            int pit1=0;
            while(p!=0) {
                PBP[pit][pit1++] = atoi(p);
                p = strtok(0, " ");
            }
            pit++;
        }
    }

	if (verbose == 1) {
		fprintf(stdout, "Forced base pairs: ");
		for(it=0; it< nFBP; it++) {
			for(int k=1;k<= FBP[it][2];k++)
				fprintf(stdout, "(%d,%d) ", FBP[it][0]+k-1, FBP[it][1]-k+1);
		}
		fprintf(stdout, "\nProhibited base pairs: ");
		for(it=0; it< nPBP; it++) {
			for(int k=1;k<= PBP[it][2];k++)
				fprintf(stdout, "(%d,%d) ", PBP[it][0]+k-1, PBP[it][1]-k+1);
		}
		fprintf(stdout, "\n\n");
	}

    return 0;
}

int init_constraints(const char* constr_file,int length) {
	load_constraints(constr_file);

	BP = (int*) malloc((length + 1) * sizeof(int));
    if (BP == NULL) {
        perror("Cannot allocate variable 'constraints'");
        exit(-1);
    }
	
	int i, it, k;

	for(i = 1; i <= length; i++) 
		BP[i] = (RNA[i]=='N')?-1:0 ;

	if (nPBP != 0) {
		for (it = 0; it < nPBP; it++) {
			if (PBP[it][2] < 1) {
				printf("Invalid entry (%d %d %d)\n", PBP[it][0], PBP[it][1],PBP[it][2]);
				continue;
			}	
			for(k= 1; k <= PBP[it][2];k++) {
				BP[PBP[it][0]+k-1] = -1;
				if(PBP[it][1]!=0) {
					BP[PBP[it][1]+1-k] = -1;
				}
			}
		}
	}
	if (nFBP != 0) {
		for (it = 0; it < nFBP; it++) {
			if (FBP[it][2] < 1) {
				printf("Invalid entry (%d %d %d)\n", FBP[it][0], FBP[it][1], FBP[it][2]);
				continue;
			}
			for(k=1; k <= FBP[it][2];k++) {
				int i1 = FBP[it][0]+k-1;
				int j1 = FBP[it][1]-k+1;
				if (!canPair(RNA[FBP[it][0]+k-1], RNA[FBP[it][1]-k+1])) {
					printf("Can't constrain (%d,%d)\n", FBP[it][0]+k-1, FBP[it][1]-k+1);
					continue;
				}
				if (j1-i1 < TURN) {
					printf("Can't constrain (%d,%d)\n", i1, j1);
					continue;
				}	
				BP[FBP[it][0]+k-1] = FBP[it][1]+1-k;
				BP[FBP[it][1]+1-k] = FBP[it][0]+k-1;
			}
		}
	}

	return 0;
}

void free_constraints(int len) {
	free(BP);
}

void print_constraints(int len) {
    int i = 1;
    for (i = 1; i <= len; i++) {
        if (BP[i] > 0 && BP[i] > i)
            printf("(");
        else if (BP[i] > 0 && BP[i] < i)
            printf(")");
        else if (BP[i] < 0)
            printf("x");
        else
            printf(".");
    }
    printf("\n");
}

int is_ss(int i, int j) {
	if (CONS_ENABLED) {
		int it;
		for (it = i + 1; it < j; it++) {
			if (BP[it] > 0) return 1;
		}
		return 0;
	}
	else
		return 0;
}

int prohibit_base(int i) {
		return (BP[i] == -1);
}

int check_base(int i) {
	if (CONS_ENABLED) 
		return (BP[i] <= 0);
	else
		return 1;
}

int force_pair(int i, int j) {
		return (BP[i] > 0 && j != BP[i]);
}

int force_pair1(int i, int j) {
	if (CONS_ENABLED) 
		return (BP[i]==j);
	else
		return 0;
}

int check_iloop(int i, int j, int p, int q) {
	if (CONS_ENABLED) 
		return is_ss(i,p) || is_ss(q,j);
	else 
		return 0;
}

int check_pair(int i, int j) {
	if (CONS_ENABLED) 
		return prohibit_base(i) || prohibit_base(j) || force_pair(i,j) || force_pair(j,i);
	else
		return 0;
}

int check_stack(int i, int j) {
	if (CONS_ENABLED) 
		return force_pair(i,j) || force_pair(j,i);
	else
		return 0;
}

int check_hairpin(int i, int j) {
	if (CONS_ENABLED) 
		return is_ss(i,j) || force_pair(i,j) || force_pair(j,i);
	else
		return 0;
}
