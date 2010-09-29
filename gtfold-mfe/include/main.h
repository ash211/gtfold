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

#ifndef _MAIN_H
#define _MAIN_H

#include <string>
using namespace std;

void traverse_wholeData(int[], int, int);

enum GTFOLD_RETURN_VAL { ERR_OPEN_FILE = 0, NO_CONS_FOUND, GTFOLD_OK};

GTFOLD_RETURN_VAL initialize_constraints(int*** fbp, int*** pbp, int& numpConstraints, int& numfConstraints, const char* constr_file);
bool handle_IUPAC_code(const std::string& s, const int bases);

void limit_contact_distance(int lCD, int length);

#endif
