#include "model.h"
#include "atom.h"
#include <iostream>
#include <fstream>
#include <string>
#include <stdlib.h>
#include <time.h>
#include <vector>
#include <tuple>
#include <map>

using namespace std;

double model::getRandomDecimal() {
    srand(time(NULL));
    return (((rand()%9999)+1)/10000.0);
}

// Constructor for initializing atoms with configAtoms first, then random atoms based off of atom_type_* and atom_number_* configuration options
model::model(map<string, tuple<int, double> > uAtoms, vector<atom> configAtoms, string efn) {
    string errorFileName = efn + ".error"; // output_file_prefix.error where output_file_prefix is from the config file
    ofstream fout;
    fout.open(errorFileName.c_str());
    if (fout.fail()) {
        cerr << "Could not open error output file " << errorFileName << endl;
    }
    
    this->atoms = configAtoms; // copy configuration atoms into the model first
    int oSize = (int)atoms.size(); // store original size of atoms in case a pseudo-clear() is necessary
    int aTotal = oSize; // Start with size of atoms, in case configAtoms was non-empty
    for(const auto &a : uAtoms) // Get total number of atoms
        aTotal += get<0>(a.second); // aTotal += quantity of each atom
    
    setBoxSize(2.715 * pow(aTotal, 1/3.0)); // sets a dimension to 2.715 * (total number of atoms)^(1/3)
    
    int numFail = 0;
    
    int id;
    while(true) {
        if(numFail >= 100) {
	    fout << "Bonding failed after 100 attempts. Program aborted." << endl;
	    if(fout)
		fout.close(); // close filestream
	    atoms.clear(); // clear atoms, as deconstructor will not be called when the exception is thrown
	    throw ("Bonding failed after 100 attempts."); // Throws an exception to stop program
	}
	
	atoms.resize(oSize); // clears randomly generated atoms
	id = 0;
	for(const auto &a : uAtoms) { // iterate through each type of atom
	    for(int i = 0; i < get<0>(a.second); i++) { // iterate through each atom of type j
		// add an atom with type "j" and unique id "id" to random coordinates in the box
		atoms.push_back(atom(getRandomDecimal()*boxSize[0], getRandomDecimal()*boxSize[1], getRandomDecimal()*boxSize[2], id, a.first));
		id++;
	    }
	}
	
	if(bondAtoms()) // If the atoms can bond, then exit the while loop
	    break;
	else
	    numFail++;
    }

    for(atom a : atoms) { // iterate through all atoms in atoms
	map<double, atom> sortedAtoms; // map to sort atoms by distance from a
	for(atom b : atoms) { // iterate through all atoms again and add to sortedAtoms
	    sortedAtoms[a.distanceTo(b)] = b;
	}
	map<double, atom>::iterator c = sortedAtoms.begin()++; // the first distance will always be 0 (distance to itself) so start at begin()++
	for(int i = 0; i < 4; i++) { // iterate through the first 4 elements for the closest 4 atoms
	    a.addNN(c->second); // add them to nn
	    c++; // iterate sortedAtoms iterator
	}
    }
    
    /*
	for(int i = 0; i < 50 && relax() > 0; i++);
	keepInBox();
    */
    fout.close();	 
}

bool model::bondAtoms() {
    for(atom a : atoms) { // iterate through all atoms in atoms
	for(atom n : a.getNN()) { // iterate through all nearest neighbors of atom a
	    if(a.getBonds() < 4 && n.getBonds() < 4) // if both atoms need more bonds, then bond
		a.bond(n);
	}
    }

    return true;
}

/* Bonds the atoms
Returns true if the bonding succeeded
bool model::bondAtoms() {
    int n = (int)atoms.size(); // get number of total atoms in the model
    
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
    //    return true;
}
*/

double model::getBoxSize() {
    return boxSize[0];
}

void model::setBoxSize(double side) {
    for(int i = 0; i < 3; i++)
	boxSize[i] = side;
}
