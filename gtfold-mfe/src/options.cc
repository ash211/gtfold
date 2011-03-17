#include "options.h"

bool ILSA;
bool NOISOLATE;
bool USERDATA;
bool PARAMS;
bool LIMIT_DISTANCE;

std::string datadir;
std::string dataparam;
std::string seqfile ;

int delta;
int nThreads;

void help() 
{
	fprintf(stderr, "\n");
	fprintf(stderr, "Usage: gtfold [OPTIONS]  FILE\n"); 
	fprintf(stderr, "\tSequence file has to be in one of the two formats: Single line or FASTA\n\n"); 
	fprintf(stderr, "OPTIONS\n" );
	fprintf(stderr, "-ni, --noisolate\n"); 
	fprintf(stderr, "\tprevents isolated base pairs from forming\n"); 
	fprintf(stderr, "-c, --constraints filename\n");
	fprintf(stderr,	"\tconstraints is a optional parameter\n");
	fprintf(stderr, "\tSyntax for giving constraints is:\n\tfor forcing (i,j)(i+1,j-1),.......,(i+k-1,j-k+1) base pair, F i j k \n");
	fprintf(stderr,	"\tto prohibit (i,j)(i+1,j-1),.......,(i+k-1,j-k+1) base pair, P i j k \n");
	fprintf(stderr, "\tto make bases from i to i+k-1 single stranded bases, P i 0 k \n");
	fprintf(stderr,	"-lcd --limitcd value\n");
	fprintf(stderr,	"\tlimits the 'contact distance' for a base pair\n");
	fprintf(stderr,	"-d --datadir dirname \n");
	fprintf(stderr, "\n");
	exit(-1);
}


void parse_options(int argc, char** argv)
{
	int fileIndex = 0;
	int consIndex = 0;
	int dataIndex = 0;
	int paramsIndex =0;
	int eIndex = 0;
	int nTIndex = 0;

	int i = 1;

	while (i < argc) {
		if (argv[i][0] == '-') {
			if (strcmp(argv[i], "-ilsa") == 0) {
				ILSA = true;
			} else if (strcmp(argv[i], "-noisolate") == 0) {
				NOISOLATE = true;
			} else if (strcmp(argv[i], "-help") == 0) {
				help();
			} else if (strcmp(argv[i], "-constraints") == 0) {
				if (i < argc)
					consIndex = ++i;
				else
					help();
			} else if (strcmp(argv[i], "-params")==0) { 
				PARAMS = true;			  
				if (i < argc)
					paramsIndex = ++i;
				else
					help();
			} else if (strcmp(argv[i], "-datadir") == 0) {
				USERDATA = true;
				if (i < argc)
					dataIndex = ++i;
				else
					help();
			} else if (strcmp(argv[i], "-e") == 0)
			{
				if (i < argc)
					eIndex = ++i;
				else
					help();	
			}
			else if (strcmp(argv[i], "-nTHREAD") == 0)
			{
				if (i < argc)
					nTIndex = ++i;
				else
					help();	
			}
		} else {
			fileIndex = i;
			seqfile = argv[fileIndex];
		}
		i++;
	}

	if (fileIndex == 0)
		help();
	
	if (NOISOLATE == true)
		fprintf(stdout, "Not allowing isolated base pairs\n");
	else
		fprintf(stdout, "Allowing isolated base pairs\n");
	
	if (consIndex != 0)
	{
		fprintf(stdout, "Constraint file index: %d\n", consIndex);
	}

	if (dataIndex != 0)
	{
		datadir = argv[dataIndex];
	}

	if (paramsIndex != 0)
	{
		dataparam = argv[paramsIndex];
	}

	nThreads = -1;
	if (nTIndex != 0)
	{
		nThreads = atoi(argv[nTIndex]);
	}
	
	delta = 0;
	if (eIndex != 0)
	{
		delta = atoi(argv[eIndex]);
	}	
	if (delta > 0)
	{
		fprintf(stdout, "Suboptimal range = %d\n\n", delta);
	}

	if (argc < 2)
	{
		help();
	}
}
