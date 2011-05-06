#ifndef _OPTIONS_H_
#define _OPTIONS_H_

#include <stdlib.h>
#include <string>
#include <cstring>
#include <stdio.h>

using namespace std;

extern bool ILSA; 
extern bool NOISOLATE;
extern bool USERDATA;
extern bool PARAMS;
extern bool LIMIT_DISTANCE;
extern bool BPP_ENABLED;
extern bool SUBOPT_ENABLED;
extern bool CONS_ENABLED;
extern bool VERBOSE;

extern string seqfile;
extern string constraintsFile;
extern string outputFile;
extern string bppOutputFile;

extern int suboptDelta;
extern int nThreads;
extern int contactDistance;

void help(); 
void parse_options(int argc, char** argv);
void printRunConfiguration(string seq);

#endif
