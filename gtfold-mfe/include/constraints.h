#ifndef _CONSTRAINTS_H_
#define _CONSTRAINTS_H_

extern int* BP;
extern int** PBP;
extern int** FBP;

extern int nPBP;
extern int nFBP;

//static int load_constraints(const char* constr_file, int verbose=0);

int init_constraints(const char* constr_file, int length) ;
void free_constraints(int length) ;
void print_constraints(int length) ;


#ifdef __cplusplus
extern "C" {
#endif
int ssOK(int i, int j); 
int baseOK(int i);
int pairOK(int i, int j);
#ifdef __cplusplus
}
#endif

#endif
