
#ifndef _ALGORITHMS_PARTITION_H
#define _ALGORITHMS_PARTITION_H

#ifdef __cplusplus
extern "C" {
#endif


void fill_partition_fn_arrays(int len, double** QB, double** Q, double** QM);
void fillBasePairProbabilities(int length, int *structure, double **Q, double **QB, double **QM, double**P);
void printBasePairProbabilities(int n, int *structure, double **P);

double probabilityUnpaired(int length, int i, double **P);

double **mallocTwoD(int r, int c);
void freeTwoD(double** arr, int r, int c);

typedef struct _pFuncData {
    int len;
    double** QB;
    double** Q;
    double** QM;

    // TODO: probability function stuff here
} pFuncData;


#ifdef __cplusplus
}
#endif

#endif
