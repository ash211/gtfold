#include "options.h"

using namespace std;

bool ILSA;
bool NOISOLATE;
bool USERDATA;
bool PARAMS;
bool LIMIT_DISTANCE;
bool ENABLE_BPP;

string seqfile = "";
string constraintsFile = "";
string outputFile = "";

int suboptDelta = -1;
int nThreads = -1;
int contactDistance = -1;

/**
 * Print the help message and quit.
 */
void help() {
    fprintf(stdout, "Usage: gtfold [OPTION]... FILE\n\n");

    fprintf(stdout, "  FILE is an RNA sequence file.  Single line or FASTA formats are accepted.\n\n");

    fprintf(stdout, "OPTIONS\n");
    fprintf(stdout, "   -c, --constraints FILE\n                        Load constraints from FILE.  See Constraint syntax below\n");
    fprintf(stdout, "   -d, --limitCD num    Set a maximum base pair contact distance to num. If no\n                        limit is given, base pairs can be over any distance\n");
    fprintf(stdout, "   -n, --noisolate      Prevent isolated base pairs from forming\n");
    fprintf(stdout, "   -o, --output FILE    Output to FILE (default output is to a .ct extension)\n");
    fprintf(stdout, "   -t, --threads num    Limit number of threads used\n");

    fprintf(stdout, "\n");
    fprintf(stdout, "   -h, --help           Output help (this message) and exit\n");
    fprintf(stdout, "   -v, --verbose        Run in verbose mode\n");

    fprintf(stdout, "\nBETA OPTIONS\n");
    fprintf(stdout, "   --bpp                Calculate base pair probabilities\n");
    fprintf(stdout, "   --subopt range       Calculate suboptimal structures within 'range' kcal/mol\n");
    fprintf(stdout, "                        of the mfe\n");

    fprintf(stdout, "\nConstraint syntax:\n\tF i j k  # force (i,j)(i+1,j-1),.......,(i+k-1,j-k+1) pairs\n\tP i j k  # prohibit (i,j)(i+1,j-1),.......,(i+k-1,j-k+1) pairs\n\tP i 0 k  # make bases from i to i+k-1 single stranded bases.\n");
    exit(-1);
}

/**
 * Parse the options from argc and argv and save them into global state.
 */
void parse_options(int argc, char** argv) {
	int i;

	for(i=1; i<argc; i++) {
		if(argv[i][0] == '-') {
			if(strcmp(argv[i], "--help") == 0 || strcmp(argv[i], "-h") == 0) {
				help();
			} else if(strcmp(argv[i], "--constraints") == 0 || strcmp(argv[i], "-c") == 0) {
				if(i < argc)
					constraintsFile = argv[++i];
				else
					help();
			} else if(strcmp(argv[i], "--limitCD") == 0 || strcmp(argv[i], "-d") == 0) {
				if(i < argc)
					contactDistance = atoi(argv[++i]);
				else
					help();
			} else if(strcmp(argv[i], "--noisolate") == 0 || strcmp(argv[i], "-n") == 0) {
				NOISOLATE = true;
			} else if(strcmp(argv[i], "--output") == 0 || strcmp(argv[i], "-o") == 0) {
				if(i < argc)
					outputFile = argv[++i];
				else
					help();
			} else if(strcmp(argv[i], "--threads") == 0 || strcmp(argv[i], "-t") == 0) {
				if(i < argc)
					nThreads = atoi(argv[++i]);
				else
					help();	
			} else if(strcmp(argv[i], "--bpp") == 0) {
				ENABLE_BPP = true;
			} else if(strcmp(argv[i], "--subopt") == 0) {
				if(i < argc)
					suboptDelta = atoi(argv[i]);
				else
					help();
			}
		} else {
			seqfile = argv[i];
		}
	}

	// Must have an input file specified
	if(seqfile.compare("") == 0) {
		help();
		printf("Missing input file.\n");
	}

	printRunConfiguration();
}

/**
 * Prints the run configuration for this run.
 *
 * The lines that start with a '-' are normal options, the '+' are beta options.
 */
void printRunConfiguration() {
	bool standardRun = true;

	printf("\nRun Configuration:\n");

#ifdef _OPENMP
	if(nThreads == -1)
		printf("- [OMP] thread count: %d\n", omp_get_num_threads());
	else
		printf("- [OMP] thread count: %d\n", nThreads);
#else
	printf("- thread count: 1\n");
#endif


	if (NOISOLATE == true) {
		printf("- preventing isolated base pairs\n");
		standardRun = false;
	}

	if(!constraintsFile.empty()) {
		printf("- using constraint file: %s\n", constraintsFile.c_str());
		standardRun = false;
	}

	if (contactDistance != -1) {
		printf("- maximum contact distance: %d\n", contactDistance);
		standardRun = false;
	}

	if (ENABLE_BPP == true) {
		printf("+ calculating base pair probabilities\n");
		standardRun = false;
	}

	if (suboptDelta != -1) {
		printf("+ calculating suboptimal structures within %d kcal/mol of MFE\n", suboptDelta);
		standardRun = false;
	}

	if(standardRun)
		printf("- standard\n");
}
