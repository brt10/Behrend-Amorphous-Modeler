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

/* Bonds the atoms
Returns true if the bonding succeeded
*/
bool model::bondAtoms()
{
	int n = atoms.size();
/*
	// this avoids dead ends where atoms in less dense areas are not bonded
	//initAtomsInShell();
	//sort(atoms.begin(), atoms.end(), isMoreDesparate);

	// select an atom
	for(int i = 0; i < n; i++)
	{
		// while it doesn't have the needed bonds, bond
		while(atoms[i].bonds.size() < NUM_BONDS[atoms[i].type]) 
		{
			size_t closestSubscript = i;
			double closestDistance = INT_MAX;
			Point closestPoint = atoms[i].coordinates;

			// all elements 0-i already have 4 bonds and can be ignored
			for(size_t potentialPartnerSubs = i+1; 
				potentialPartnerSubs < n; 
				potentialPartnerSubs++)
			{

				// a partner is valid if
				// 1. this partner does not have 4 bonds
				// 2. this partner is not already bonded to our current atom
				if(atoms[potentialPartnerSubs].bonds.size() < NUM_BONDS[atoms[potentialPartnerSubs].type] &&
					find(atoms[i].bonds.begin(), atoms[i].bonds.end(), 
					Partner(Point(0,0,0), atoms[potentialPartnerSubs].id)) == atoms[i].bonds.end())
				{
					// find closest version of this potential partner
					Point p = closestPointWithPBC(atoms[i].coordinates, atoms[potentialPartnerSubs].coordinates);
					double tempDistance = p.distanceTo(atoms[i].coordinates);

					// if it is the closest possible partner, take note
					if(tempDistance < closestDistance)
					{
						closestDistance = tempDistance;
						closestSubscript = potentialPartnerSubs;
						closestPoint = p;
					}
				}
			}

			// Theoretically, this should not happen
			if(closestSubscript == i)
				return false;

			// adjust bonds
			atoms[i].bonds.push_back(Partner(closestPoint, atoms[closestSubscript].id, closestDistance));

			// this takes care of a periodic boundary conditions problem
			Vector displacement(closestPoint.x - atoms[i].coordinates.x,
				closestPoint.y - atoms[i].coordinates.y,
				closestPoint.z - atoms[i].coordinates.z);

			atoms[closestSubscript].bonds.push_back(Partner(Point(atoms[closestSubscript].coordinates.x - displacement.x,
				atoms[closestSubscript].coordinates.y - displacement.y,
				atoms[closestSubscript].coordinates.z - displacement.z),
				atoms[i].id, closestDistance));
		}
	}

	// resorts the atoms by ID so that atom #n will be at index n-1 in the array
	sort(atoms.begin(), atoms.end(), hasSmallerID);
*/
	return true;
}
