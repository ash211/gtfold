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

#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <math.h>
#include <time.h>
#include <sys/time.h>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "main.h"
#include "utils.h"
#include "loader.h"
#include "options.h"
#include "global.h"
#include "energy.h"
#include "algorithms.h"
#include "traceback.h"
#include "subopt_traceback.h"

using namespace std;

double get_seconds() 
{
	struct timeval tv;
	gettimeofday(&tv, NULL);
	return (double) tv.tv_sec + (double) tv.tv_usec / 1000000.0;
}

void init_fold(int len) 
{	

 	init_global_params(len);
	create_tables(len);
	return;
}

void free_fold(int len) 
{
	free_global_params();
	free_tables(len);
	return;
}

/**
 * Read the sequence out of the given filename and store it in seq
 *
 * @param filename A c string with the file to open
 * @param seq A C++ string object (passed by reference) to save to
 * @return SUCCESS or FAILURE
 */
int read_sequence_file(const char* filename, std::string& seq)
{
	seq = "";

	ifstream fs;
	fs.open(filename, ios::in);
	if (fs == NULL)
		return FAILURE;

	string line;
	getline(fs, line);
	while(line.length() > 0) {
		// exclude lines starting with FASTA comment characters
		if(line[0] != ';' && line[0] != '>')
			seq += line;
		getline(fs, line);
	}

	fs.close();

	return SUCCESS;
}

int handle_IUPAC_code(const char* str, const int bases)
{
	int* stack_unidentified_base;
	int stack_count=0;
	bool unspecd=0;
	stack_unidentified_base=new int[bases];
	std::string s = str;

	for(int i = 1; i <= bases; i++) 
	{
		RNA[i] = getBase((s.substr(i-1,1)).c_str());
		RNA1[i] = getBase1((s.substr(i-1,1)).c_str());

		if (RNA[i]=='X') 
		{
			fprintf(stderr,"ERROR: Base unrecognized\n");
			return FAILURE;
		}
		else if(RNA[i]!='X' && RNA1[i]=='N')
		{
			unspecd=1;
			stack_unidentified_base[stack_count]=i;
			stack_count++;
		}
	}
	if(unspecd) 
	{
		printf("IUPAC codes have been detected at positions:");

		for(int i=0;i<stack_count;i++)
		{
			printf("%d , ",stack_unidentified_base[i]);
		}
		printf("\n");
		//printf("You may wish to resubmit the sequence with fully specified positions. Alternatively, GTfold will fold the sequence under the standard assumption that these ambiguous positions do not pair.  Do you wish to continue with the current computation? <Y/N>");
		//char reply;
		//scanf("%c",&reply);
		//return (reply=='n'||reply=='N')?(FAILURE):(SUCCESS);
		return SUCCESS;
	}
	else 
	{
		return SUCCESS;
	}
}

int main(int argc, char** argv) {
	std::string seq;
	int energy;
	double t1;

	fprintf(stdout,"GTfold: A Scalable Multicore Code for RNA Secondary Structure Prediction\n");
	fprintf(stdout,"(c) 2007-2010  D.A. Bader, S. Mallidi, A. Mathuriya, C.E. Heitsch, S.C. Harvey\n");
	fprintf(stdout, "Georgia Institute of Technology\n");
	
	parse_options(argc, argv);

	fprintf(stdout, "Opening file: %s\n", seqfile.c_str());
	if (read_sequence_file(seqfile.c_str(), seq) == 0)
	{
		fprintf(stdout, "File open failed.\n\n");
		exit(-1);
	}
	fprintf(stdout, "Sequence length: %5d\n", seq.length());
	fprintf(stdout, "Sequence : %s\n", seq.c_str());
	
	init_fold(seq.length());
	
	if (handle_IUPAC_code(seq.c_str(), seq.length())  == 0)
	{
		free_fold(seq.length());
		exit(0);
	}
		
	if(USERDATA==true)
		populate(datadir.c_str(),true);
	else if (PARAMS == true)
		populate(dataparam.c_str(),false);
	else
		populate("Turner99",false); 
	
	fprintf(stdout,"Computing minimum free energy structure. . . \n");
	fflush(stdout);

	t1 = get_seconds();
	energy = calculate(seq.length());
	t1 = get_seconds() - t1;
	
	fprintf(stdout,"Done.\n\n");
	fprintf(stdout,"Minimum Free Energy = %12.4f\n", energy/100.00);
	fprintf(stdout,"MFE running time (in seconds): %9.6f\n\n", t1);
	
	
	if (delta > 0)
	{	
		t1 = get_seconds();
		subopt_traceback(seq.length(), delta);
		t1 = get_seconds() - t1;
		fprintf(stdout,"Traceback running time (in seconds): %9.6f\n\n", t1);

		free_fold(seq.length());
		exit(0);
	}
	
	t1 = get_seconds();
	trace(seq.length());
	t1 = get_seconds() - t1;
	
	/*
	std::stringstream ss1, ss2;
	char suboptfile[1024];
	ss1 << seq.length();
	ss2 << energy/100.0;

	strcpy (suboptfile, argv[fileIndex]);
	strcat ( suboptfile, ".ct");
	ofstream outfile;
	outfile.open ( suboptfile );

	fprintf(stdout, "Writing secondary structure to the file: %s\n", suboptfile);

	outfile << bases << "\t  dG = " << energy/100.0;
	i = 1;
	while ( i <= bases ) {
		outfile << endl << i << "\t" << s[i-1] << "\t" << i-1 << "\t" << (i+1)%(bases+1) << "\t" << structure[i] << "\t" << i;
		i++;
	}
	outfile << endl;

	outfile.close();

	fprintf(stdout,"\n\n");
	fprintf(stdout,"Traceback running time (in seconds): %9.6f\n", t1);

	fprintf(stdout, "\n\nFolding complete\n\n");
	printSequence(bases);
	//printConstraints(bases);
	printStructure(bases);
*/
	free_fold(seq.length());

	return 0;
}


int initialize_constraints(int*** fbp, int ***pbp, int& numpConstraints, int& numfConstraints, const char* constr_file)
{
	ifstream cfcons;

	fprintf(stdout, "Running with constraints\n");
	fprintf(stdout, "Opening constraint file: %s\n", constr_file);

	cfcons.open(constr_file, ios::in);
	if (cfcons != NULL)
		fprintf(stdout, "Constraint file opened.\n");
	else {
		fprintf(stderr, "Error opening constraint file\n\n");
		cfcons.close();
		return FAILURE; //exit(-1);
	}

	char cons[100];

	while (!cfcons.eof()) {
		cfcons.getline(cons, 100);
		if (cons[0] == 'F' || cons[0] == 'f')
			numfConstraints++;
		if (cons[0] == 'P' || cons[0] == 'p')
			numpConstraints++;
	}
	cfcons.close();

	fprintf(stdout, "Number of Constraints given: %d\n\n", numfConstraints
			+ numpConstraints);
	if (numfConstraints + numpConstraints != 0)
		fprintf(stdout, "Reading Constraints.\n");
	else {
		fprintf(stderr, "No Constraints found.\n\n");
		return FAILURE;
	}

	*fbp = (int**) malloc(numfConstraints * sizeof(int*));
	*pbp = (int**) malloc(numpConstraints * sizeof(int*));

	int fit = 0, pit = 0, it = 0;

	for (it = 0; it < numfConstraints; it++) {
		(*fbp)[it] = (int*) malloc(2* sizeof (int));
	}
	for(it=0; it<numpConstraints; it++) {
		(*pbp)[it] = (int*)malloc(2*sizeof(int));
	}
	cfcons.open(constr_file, ios::in);

	while(!cfcons.eof()) {
		cfcons.getline(cons,100);
		char *p=strtok(cons, " ");
		p = strtok(NULL, " ");
		if(cons[0]=='F' || cons[0]=='f') {
			int fit1=0;
			while(p!=NULL) {
				(*fbp)[fit][fit1++] = atoi(p);
				p = strtok(NULL, " ");
			}
			fit++;
		}
		if( cons[0]=='P' || cons[0]=='p') {
			int pit1=0;
			while(p!=NULL) {
				(*pbp)[pit][pit1++] = atoi(p);
				p = strtok(NULL, " ");
			}
			pit++;
		}
	}

	fprintf(stdout, "Forced base pairs: ");
	
	for(it=0; it<numfConstraints; it++) 
	{
		for(int k=1;k<= (*fbp)[it][2];k++)
			fprintf(stdout, "(%d,%d) ", (*fbp)[it][0]+k-1, (*fbp)[it][1]-k+1);
	}
	fprintf(stdout, "\nProhibited base pairs: ");
	for(it=0; it<numpConstraints; it++) 
	{
		for(int k=1;k<= (*pbp)[it][2];k++)
			fprintf(stdout, "(%d,%d) ", (*pbp)[it][0]+k-1, (*pbp)[it][1]-k+1);
	}

	fprintf(stdout, "\n\n");

	return SUCCESS;
}
