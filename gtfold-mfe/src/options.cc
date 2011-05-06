#include "loader.h"
#include "options.h"

using namespace std;

bool ILSA;
bool NOISOLATE;
bool USERDATA;
bool PARAMS;
bool LIMIT_DISTANCE;
bool BPP_ENABLED;
bool SUBOPT_ENABLED;
bool CONS_ENABLED = false;
bool VERBOSE = false;

string seqfile = "";
string constraintsFile = "";

string outputFileBase = "";
string outputFile = "";
string bppOutputFile = "";

int suboptDelta = -1;
int nThreads = -1;
int contactDistance = -1;

/**
 * Print the help message and quit.
 */
void help() {
    printf("Usage: gtfold [OPTION]... FILE\n\n");

    printf("  FILE is an RNA sequence file.  Single line or FASTA formats are accepted.\n\n");

    printf("OPTIONS\n");
    printf("   -c, --constraints FILE\n");
    printf("                        Load constraints from FILE.  See Constraint syntax below\n");
    printf("   -d, --limitCD num    Set a maximum base pair contact distance to num. If no\n");
    printf("                        limit is given, base pairs can be over any distance\n");
    printf("   -n, --noisolate      Prevent isolated base pairs from forming\n");
    printf("   -o, --output FILE    Output to FILE (default output is to a .ct extension)\n");
    printf("   -t, --threads num    Limit number of threads used\n");

    printf("\n");
    printf("   -h, --help           Output help (this message) and exit\n");
    printf("   -v, --verbose        Run in verbose mode\n");

    printf("\nBETA OPTIONS\n");
    printf("   --bpp                Calculate base pair probabilities\n");
    printf("   --subopt range       Calculate suboptimal structures within 'range' kcal/mol\n");
    printf("                        of the mfe\n");

    printf("\nConstraint syntax:\n");
    printf("\tF i j k  # force (i,j)(i+1,j-1),.......,(i+k-1,j-k+1) pairs\n");
    printf("\tP i j k  # prohibit (i,j)(i+1,j-1),.......,(i+k-1,j-k+1) pairs\n");
    printf("\tP i 0 k  # make bases from i to i+k-1 single stranded bases.\n");
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
				if(i < argc) {
					constraintsFile = argv[++i];
					CONS_ENABLED = true;
				}
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
			} else if (strcmp(argv[i], "--verbose") == 0 || strcmp(argv[i], "-v") == 0) {
                VERBOSE = true;
			}
			else if(strcmp(argv[i], "--bpp") == 0) {
				BPP_ENABLED = true;
			} else if(strcmp(argv[i], "--subopt") == 0) {
				SUBOPT_ENABLED = true;
				if(i < argc)
					suboptDelta = atoi(argv[++i]);
				else
					help();
			}
		} else {
			seqfile = argv[i];
		}
	}

	// Must have an input file specified
	if(seqfile.empty()) {
		help();
		printf("Missing input file.\n");
	}

    // Create the extension-less outputFileBase:

    // prefer to base it on the specified output file
    outputFileBase = "";
    if(!outputFile.empty())
        outputFileBase += outputFile.c_str();
    else
        outputFileBase += seqfile.c_str();

    // drop any existing extension
    if(outputFileBase.find(".") != string::npos)
        outputFileBase.erase(outputFileBase.rfind("."));

    // append .ct to get regular outputFile
    outputFile = outputFileBase + ".ct";

    // append .bpp to get bppOutputFile
    bppOutputFile = outputFileBase + ".bpp";
}

/**
 * Prints the run configuration for this run.
 *
 * The lines that start with a '-' are normal options, the '+' are beta options.
 */
void printRunConfiguration(string seq) {
	bool standardRun = true;

	printf("Run Configuration:\n");

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

	if (BPP_ENABLED == true) {
		printf("+ calculating base pair probabilities\n");
		standardRun = false;
	}

	if (SUBOPT_ENABLED) {
		printf("+ calculating suboptimal structures within %d kcal/mol of MFE\n", suboptDelta);
		standardRun = false;
	}

	if(standardRun)
		printf("- standard\n");

	printf("- thermodynamic parameters: %s\n", EN_DATADIR.c_str());
	printf("- input file: %s\n", seqfile.c_str());
	printf("  - sequence length: %d\n", (int)seq.length());
	printf("  - sequence contents: %s\n", seq.c_str());
	printf("- output file: %s\n", outputFile.c_str());
}
