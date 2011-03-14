#ifndef _OPTIONS_H_
#define _OPTIONS_H_

#include <stdlib.h>
#include <string>
#include <cstring>
#include <stdio.h>

extern bool ILSA; 
extern bool NOISOLATE;
extern bool USERDATA;
extern bool PARAMS;
extern bool LIMIT_DISTANCE;

extern std::string datadir;
extern std::string dataparam;
extern std::string seqfile ;

extern int delta;
extern int nThreads;

void help(); 
void parse_options(int argc, char** argv);

#endif
