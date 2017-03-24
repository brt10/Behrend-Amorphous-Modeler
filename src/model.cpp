#include "model.h"
#include "atom.h"
#include <iostream>
#include <fstream>
#include <string>

using namespace std;

model::model() {
	// placeholder for model default constructor
}

double model::getRandomDecimal() {
	return (((rand()%9999)+1)/10000.0);
}

bool model::makeDefectedSystem() {
	string errorFileName = "output.error"; //FILE_NAME + ".error";  string for the error file name
	ofstream fout;
	fout.open(errorFileName.c_str());
	if (fout.fail()) {
		cerr << "Could not open error output file " << errorFileName << endl;
		return 1;
	}

	boxSize[0] = 1000.0; //2.715 * pow(NEW_ATOMS, 1/3.0); // sets a dimension to 2.715 * NEW_ATOMS^(1/3) ?
	boxSize[2] = boxSize[1] = boxSize[0]; // a box

	int numFail = 0;

	while(true) {
		if(numFail >= 100){
			fout << "Bonding failed after 100 attempts. Program aborted." << endl;
			return 1;
		}

		atoms.clear(); // clears atoms vector
		int id = 1; // unique integer id for each atom, probably redundant because vector indices are unique

		for(int j = 0; j < 3 /*N_TYPES*/; j++){ // iterate through each type of atom
			for(int i = 0; i < 10 /*nAtoms[j]*/; i++) { // iterate through each atom of type j
				atoms.push_back(atom(getRandomDecimal()*boxSize[0], getRandomDecimal()*boxSize[1], getRandomDecimal()*boxSize[2], id, 0 /*nTypes[j]*/));
				id++;
			}
		}

		// Can't find a bondAtoms() function anywhere in the project
		/*
		if(bondAtoms()){
			break;
		}
		else{
			numFail++;
		}
		*/
	}

	// not implemented yet. initClosestNeighbors clears vector<nn> in each atom of atoms[] then finds the correct nn
/*	initClosestNeighbors();
	for(int i = 0; i < 50 && relax() > 0; i++);
	keepInBox();
*/
	fout.close();

	return 0;
}
