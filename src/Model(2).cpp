/* Alicia Klinvex
April 22, 2008

Last modified by Brittany Hoard on June 16, 2012

Last modified by BRT in 2026 using co-pilot

Stores a set of atoms
Contains functions for amorphizing the atoms
*/
//#include "mpi.h"
#include "Model.h"
//#include "pthread.h"
//#include "Python.h"
#include <algorithm>
#include <fstream>
#include <iostream>
#include <ctime>
#include <queue>
#include <cassert>
#include <iomanip>
#include <string>
#include <sstream>
#include <math.h>
#include <climits>

/* Returns true if atom a is in a sparser area than atom b
*/
bool isMoreDesparate(const Atom &a, const Atom &b)
{
	return (a.atomsInShell < b.atomsInShell);
}

/* Returns true if atom a's ID is less than atom b's ID
*/
bool hasSmallerID(const Atom &a, const Atom &b)
{
	return (a.id < b.id);
}

/* Returns the minimum of a and b
*/
double min(double a, double b)
{
	if(a < b)
		return a;
	return b;
}

/* Returns the maximum of a and b
*/
double max(double a, double b)
{
	if(a > b)
		return a;
	return b;
}

/* Default constructor
*/
Model::Model()
{
	// seeds the random number generator
	srand((unsigned int)time(0));
	atoms.clear();

	IGNORE_ATOMS = 0;  //ignore the first n atoms from swaps
	FIX_ATOMS = 0;     //fix the first n atoms' position
	P_BONDS = 0;       //print bond data
	P_VASP = 0;        //print the model in .vasp format
	WRITE_NN = 0;      //print the closest neighbors for every atom
	N_TYPES = 0;
	NEW_ATOMS = 0;
	CREATE_MODEL = 0;  //create a new model and output to a new file
	RELAX_V = 0;
	RELAX_V_SWITCHES = 0;
	SI_3 = 0;
	SI_4 = 0;
	SI_5 = 0;
	NUM_H = 0;
	
	SHELL_ONLY = false;
	INTER = false;
	BOX_SIZE[0] = 10.86;
	BOX_SIZE[1] = 10.86;
	BOX_SIZE[2] = 10.86;
	OLD_X = 0.0;

	SHELL_SIZE = 4;

	framePause = 5;
	kT = 1;
	I_NUM_SWITCHES = 0;
	I_RELAX_V_SWITCHES = 0;
	REPEL_DISTANCE = 7.2;
	FORCE_CONVERSION = 1.602E-29;
	CLUSTER_MAX_SHELL = 4;

	RADIAL_ON = true;
	ANGULAR_ON = true;
	REPEL_ON = true;

	setDefaultConstants();
}
//
//Model::~Model()
//{
//	if(Py_IsInitialized())
//
//}

/* Bonds the atoms
Returns true if the bonding succeeded
*/
bool Model::bondAtoms()
{
	size_t n = atoms.size();

	// this avoids dead ends where atoms in less dense areas are not bonded
	initAtomsInShell();
	sort(atoms.begin(), atoms.end(), isMoreDesparate);

	// select an atom
	for(size_t i = 0; i < n; i++)
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
	
	return true;
}

/*Change the box size (for a cubic box) of an existing model, preserving the atoms' fractional coordinates
*/
void Model::changeBoxSize(double newBoxSize){

	//Change the atoms' x, y, and z coordinates
	for(size_t i = 0; i < atoms.size(); i++)
	{
		atoms[i].coordinates.x = newBoxSize*(atoms[i].coordinates.x/BOX_SIZE[0]);
		atoms[i].coordinates.y = newBoxSize*(atoms[i].coordinates.y/BOX_SIZE[1]);
		atoms[i].coordinates.z = newBoxSize*(atoms[i].coordinates.z/BOX_SIZE[2]);

		for(size_t j = 0; j < atoms[i].bonds.size(); j++){
			atoms[i].bonds[j].coordinates.x = newBoxSize*(atoms[i].bonds[j].coordinates.x/BOX_SIZE[0]);
			atoms[i].bonds[j].coordinates.y = newBoxSize*(atoms[i].bonds[j].coordinates.y/BOX_SIZE[1]);
		    atoms[i].bonds[j].coordinates.z = newBoxSize*(atoms[i].bonds[j].coordinates.z/BOX_SIZE[2]);
		}

	}
	
	//Change box size
	BOX_SIZE[0] = newBoxSize;
	BOX_SIZE[1] = newBoxSize;
	BOX_SIZE[2] = newBoxSize;

	//Relax atoms
	initClosestNeighbors();
	for(int i = 0; i < 50 && relax() > 0; i++);
	if(INTER==true){
			keepInInterfaceBox(OLD_X);
	}
	else{
			keepInBox();
	}

	return;
}		

/* Finds the version of atom p2 that's closest to atom p1
*/
Point Model::closestPointWithPBC(const Point &p1, const Point &p2)
{

	static Point bestPoint;
	/*int dx, dy, dz;
	int startX = int(p2.x / BOX_SIZE[0]) - 1, 
		startY = int(p2.y / BOX_SIZE[1]) - 1, 
		startZ = int(p2.z / BOX_SIZE[2]) - 1;
	double minDistance;
	Point temp, bestPoint;

	if(p2.x < 0)
		startX--;
	if(p2.y < 0)
		startY--;
	if(p2.z < 0)
		startZ--;

	for(dx = startX, minDistance=INT_MAX; dx <= 1; dx++)
	{
		temp.x = p2.x + dx * BOX_SIZE[0]; //x box size

		if(fabs(temp.x - p1.x) < minDistance)
		{
			minDistance = fabs(temp.x - p1.x);
			bestPoint.x = temp.x;
		}
	}

	for(dy = startY, minDistance=INT_MAX; dy <= 1; dy++)
	{
		temp.y = p2.y + dy * BOX_SIZE[1];//y box size

		if(fabs(temp.y - p1.y) < minDistance)
		{
			minDistance = fabs(temp.y - p1.y);
			bestPoint.y = temp.y;
		}
	}

	for(dz = startZ, minDistance=INT_MAX; dz <= 1; dz++)
	{
		temp.z = p2.z + dz * BOX_SIZE[2];//z box size

		if(fabs(temp.z - p1.z) < minDistance)
		{
			minDistance = fabs(temp.z - p1.z);
			bestPoint.z = temp.z;
		}
	}
	*/
	bestPoint.x = p2.x - round( (p2.x-p1.x) / BOX_SIZE[0] )*BOX_SIZE[0];
	bestPoint.y = p2.y - round( (p2.y-p1.y) / BOX_SIZE[1] )*BOX_SIZE[1];
	bestPoint.z = p2.z - round( (p2.z-p1.z) / BOX_SIZE[2] )*BOX_SIZE[2];
	return bestPoint;
}

/* computes the cosine of an angle using the dot product
*/
double Model::cosijk(const Point &i, const Point &j, const Point &k,
					 const double rij, const double rik)
{
	return (((j.x - i.x) * (k.x - i.x) + 
		(j.y - i.y) * (k.y - i.y) + 
		(j.z - i.z) * (k.z - i.z)) / 
		(rij * rik));
}

/* Sets and returns the energy of the system
*/
double Model::getEnergy()
{
	initClosestNeighbors();

	double energy = 0;
	if(RADIAL_ON)	setRadialEnergy();
	if(ANGULAR_ON)	setAngularEnergy();
	if(REPEL_ON)	setRepulsiveEnergy();

	double rEnergy = 0;
	double rpEnergy = 0;
	double angEnergy = 0;

	for(size_t i = 0; i < atoms.size(); i++){
		energy = energy +  atoms[i].radialEnergy + atoms[i].angularEnergy + atoms[i].repulsiveEnergy;
		rEnergy = rEnergy + atoms[i].radialEnergy;
		rpEnergy = rpEnergy + atoms[i].repulsiveEnergy;
		angEnergy = angEnergy + atoms[i].angularEnergy;
	}
	//cout << "Radial + Angular = " << (rEnergy+angEnergy) << endl;
	return energy;
}

/* Prints initial energy of the system
	Outputs in format
	atomID rad_e ang_e tot_e
*/
void Model::writeInitialEnergy(string file_name)
{
	initClosestNeighbors();
	if(RADIAL_ON)	setRadialEnergy();
	if(ANGULAR_ON)	setAngularEnergy();
	if(REPEL_ON)	setRepulsiveEnergy();

	ofstream eout;
	eout.open(file_name.c_str());	

	double rEnergy = 0;
	double rpEnergy = 0;
	double angEnergy = 0;

	for(size_t i = 0; i < atoms.size(); i++){		
		eout << i+1 << " " << atoms[i].radialEnergy << " "
			<< atoms[i].angularEnergy << " "
			<< (atoms[i].radialEnergy+atoms[i].angularEnergy)
			<< endl;
	}
}

/* Prints initial forces of the system
	Outputs in format
	atomID f_x f_y f_z
*/
void Model::writeInitialForces(string file_name)
{
	// SLT - Struct used to separate the radial and angular components for viewing
	struct components{
	double x, y, z;
	};
	struct forcePrint{
	components radial,angular;
	};



	size_t n = atoms.size();
	Partner temp;

	size_t i, j, k, jSubs, kSubs;
	double multiplier, tempNum, cosijk, tempJ, tempK, dist;

	Point p;
	// SLT
	//forcePrint output[8];
        vector<forcePrint> output(atoms.size());

	// reset all forces
	for(i = 0; i < n; i++)
	{
		if((SHELL_ONLY && atoms[i].shellNumber > CLUSTER_MAX_SHELL))	
			//this line was changed as a test
			continue;

		output[i].radial.x = 0; output[i].angular.x = 0;
		output[i].radial.y = 0; output[i].angular.y = 0;
		output[i].radial.z = 0; output[i].angular.z = 0;
	}

	for(i = 0; i < n; i++)
	{
		if((SHELL_ONLY && atoms[i].shellNumber > CLUSTER_MAX_SHELL))	
			continue;

		for(j = 0; j < atoms[i].bonds.size(); j++)
		{
			jSubs = atoms[i].bonds[j].id-1;

			// set radial forces of two partner atoms at a time
			// twice as efficient, since it takes the "equal and opposite" rule into account
			if(RADIAL_ON  &&  i < jSubs)
			{
				multiplier = (K_B[atoms[i].type][atoms[atoms[i].bonds[j].id-1].type] * (atoms[i].bonds[j].dist-B_0[atoms[i].type][atoms[atoms[i].bonds[j].id-1].type]) / atoms[i].bonds[j].dist);

				tempNum = multiplier * (atoms[i].bonds[j].coordinates.x - atoms[i].coordinates.x);
				output[i].radial.x += tempNum;
				output[jSubs].radial.x -= tempNum;

				tempNum = multiplier * (atoms[i].bonds[j].coordinates.y - atoms[i].coordinates.y);
				output[i].radial.y += tempNum;
				output[jSubs].radial.y -= tempNum;

				tempNum = multiplier * (atoms[i].bonds[j].coordinates.z - atoms[i].coordinates.z);
				output[i].radial.z += tempNum;
				output[jSubs].radial.z -= tempNum;				

			} // end if

			// set angular force of three atoms at a time
			// i is the "middle" atom and j/k are endpoints of the angle
			for(k = j+1; k < atoms[i].bonds.size() && ANGULAR_ON; k++)
			{
				kSubs = atoms[i].bonds[k].id-1;

				// avoids recomputing these repeatedly
				cosijk = ((atoms[i].bonds[j].coordinates.x - atoms[i].coordinates.x) * (atoms[i].bonds[k].coordinates.x - atoms[i].coordinates.x) + 
					(atoms[i].bonds[j].coordinates.y - atoms[i].coordinates.y) * (atoms[i].bonds[k].coordinates.y - atoms[i].coordinates.y) + 
					(atoms[i].bonds[j].coordinates.z - atoms[i].coordinates.z) * (atoms[i].bonds[k].coordinates.z - atoms[i].coordinates.z)) / 
					(atoms[i].bonds[j].dist * atoms[i].bonds[k].dist);
				multiplier = K_OMEGA[atoms[i].type][atoms[atoms[i].bonds[j].id-1].type][atoms[atoms[i].bonds[k].id-1].type] * 
					(COS_ANGLE[atoms[i].type]-cosijk);

				// find force on j
				tempJ = multiplier * ((atoms[i].bonds[k].coordinates.x-atoms[i].coordinates.x) / atoms[i].bonds[k].dist +
					(atoms[i].coordinates.x-atoms[i].bonds[j].coordinates.x) * cosijk / atoms[i].bonds[j].dist) / atoms[i].bonds[j].dist;
				output[jSubs].angular.x += tempJ;
				output[i].angular.x -= tempJ;

				tempJ = multiplier * ((atoms[i].bonds[k].coordinates.y-atoms[i].coordinates.y) / atoms[i].bonds[k].dist +
					(atoms[i].coordinates.y-atoms[i].bonds[j].coordinates.y) * cosijk / atoms[i].bonds[j].dist) / atoms[i].bonds[j].dist;
				output[jSubs].angular.y += tempJ;
				output[i].angular.y -= tempJ;

				tempJ = multiplier * ((atoms[i].bonds[k].coordinates.z-atoms[i].coordinates.z) / atoms[i].bonds[k].dist +
					(atoms[i].coordinates.z-atoms[i].bonds[j].coordinates.z) * cosijk / atoms[i].bonds[j].dist) / atoms[i].bonds[j].dist;
				output[jSubs].angular.z += tempJ;
				output[i].angular.z -= tempJ;

				// find force on k
				tempK = multiplier * ((atoms[i].bonds[j].coordinates.x-atoms[i].coordinates.x) / atoms[i].bonds[j].dist +
					(atoms[i].coordinates.x-atoms[i].bonds[k].coordinates.x) * cosijk / atoms[i].bonds[k].dist) / atoms[i].bonds[k].dist;
				output[kSubs].angular.x += tempK;
				output[i].angular.x -= tempK;

				tempK = multiplier * ((atoms[i].bonds[j].coordinates.y-atoms[i].coordinates.y) / atoms[i].bonds[j].dist +
					(atoms[i].coordinates.y-atoms[i].bonds[k].coordinates.y) * cosijk / atoms[i].bonds[k].dist) / atoms[i].bonds[k].dist;
				output[kSubs].angular.y += tempK;
				output[i].angular.y -= tempK;

				tempK = multiplier * ((atoms[i].bonds[j].coordinates.z-atoms[i].coordinates.z) / atoms[i].bonds[j].dist +
					(atoms[i].coordinates.z-atoms[i].bonds[k].coordinates.z) * cosijk / atoms[i].bonds[k].dist) / atoms[i].bonds[k].dist;
				output[kSubs].angular.z += tempK;
				output[i].angular.z -= tempK;
			} // end k loop
		} // end j loop

	} // end i loop


	
	ofstream fout;
	fout.open(file_name.c_str());

	//setForces();	//Compute forces, initially

	for(int i=0;i<atoms.size();i++){
		fout << i+1 << endl 
			<< " Radial " << output[i].radial.x << " "
			<< output[i].radial.y << " "
			<< output[i].radial.z << endl
			<< " Angular " << output[i].angular.x << " "
			<< output[i].angular.y << " "
			<< output[i].angular.z << endl;
	}
}

/* Returns the largest force on an atom in the system
*/
double Model::getLargestForce()
{
	size_t n = atoms.size();
	double maxForce = 0;

	for(size_t i = 0; i < n; i++)
	{
		if((SHELL_ONLY && atoms[i].shellNumber > CLUSTER_MAX_SHELL))	
			continue;

		// using absolute value
		maxForce = max(maxForce, fabs(atoms[i].force.x));
		maxForce = max(maxForce, fabs(atoms[i].force.y));
		maxForce = max(maxForce, fabs(atoms[i].force.z));
	}

	return maxForce;
}

/* Returns the metropolis acceptance probability
*/
double Model::getMetropolisProb(const double Einitial,
                                const double Efinal,
                                const double kT_local)
{
    return min(1, exp((Einitial - Efinal)/kT_local));
}
/* Returns a random decimal in the range (0, 1)
*/
double Model::getRandomDecimal()
{
	return (((rand()%9999)+1)/10000.0);
}

/* Returns a random integer in the range (1, 100)
*/
unsigned int Model::getRandomInteger()
{
	return ((rand()%100)+1);
}

/* Finds the number of atoms in each atom's shell
i.e. the sparsity of the area surrounding a particular atom
*/
void Model::initAtomsInShell()
{
	size_t n = atoms.size();

	for(size_t i = 0; i < n; i++)
	{
		for(size_t j = 0; j < n; j++)
		{
			Point p = closestPointWithPBC(atoms[i].coordinates, atoms[j].coordinates);

			if(p.distanceTo(atoms[i].coordinates) <= SHELL_SIZE)
				atoms[i].atomsInShell++;
		}
	}
}

/* Initializes the atoms nearby contributing to the repulsive energy
*/
void Model::initClosestNeighbors()
{
	double tempDist;

	for(size_t i = 0; i < atoms.size(); i++)
	{
		if(SHELL_ONLY && atoms[i].shellNumber > CLUSTER_MAX_SHELL+1)	continue;

		// reset array of nearby atoms
		atoms[i].closestNeighbors.clear();

		for(size_t j = 0; j < atoms.size(); j++)
		{
			if(SHELL_ONLY && atoms[j].shellNumber > CLUSTER_MAX_SHELL+1)	continue;

			// ignore the atom if it is a neighbor of the current one
			if(i == j || isNearestNeighbor(j+1, atoms[i]) || isSecondNearestNeighbor(j+1, i))
				continue;

			tempDist = atoms[i].coordinates.distanceTo(closestPointWithPBC(atoms[i].coordinates, atoms[j].coordinates));

			if(tempDist <= REPEL_DISTANCE)
				atoms[i].closestNeighbors.push_back(j+1);
		}
	}
}

/* Determines how far from a bond transposition all the atoms are
Uses a modified breadth first search
*/
void Model::initShellNumbers(const size_t aSubs, const size_t bSubs, 
							 const size_t cSubs, const size_t dSubs)
{
	// reset shell numbers
	for(size_t i = 0; i < atoms.size(); i++)
		atoms[i].shellNumber = CLUSTER_MAX_SHELL + 10;

	// set initial shell numbers of transposition atoms
	atoms[aSubs].shellNumber = 0;
	atoms[bSubs].shellNumber = 0;
	atoms[cSubs].shellNumber = 0;
	atoms[dSubs].shellNumber = 0;

	queue<size_t> ids;

	ids.push(aSubs+1);
	ids.push(bSubs+1);
	ids.push(cSubs+1);
	ids.push(dSubs+1);

	while(!ids.empty())
	{
		// remove the first element
		size_t tempID = ids.front();
		ids.pop();

		// update partners if current element is not on edge of cluster
		if(atoms[tempID-1].shellNumber <= CLUSTER_MAX_SHELL)
		{
			for(size_t i = 0; i < atoms[tempID-1].bonds.size(); i++)
			{
				int tempBondID = atoms[tempID-1].bonds[i].id;

				// if it hasn't been tested, push it in the queue
				if(atoms[tempBondID-1].shellNumber > CLUSTER_MAX_SHELL)
					ids.push(tempBondID);

				// calculate the shell number
				atoms[tempBondID-1].shellNumber = min(atoms[tempBondID-1].shellNumber,
					atoms[tempID-1].shellNumber+1);
			}
		}
	}
}

/* Returns true if atoms[id-1] is bonded to atom
*/
bool Model::isNearestNeighbor(const size_t id, const Atom &atom)
{
	for(size_t i = 0; i < atom.bonds.size(); i++)
	{
		if(atom.bonds[i].id == id)
			return true;
	}
	return false;
}

/* Returns true if there exists an atom p such that
atoms[id-1] is bonded to p and atoms[currentSubs] is bonded to p
*/
bool Model::isSecondNearestNeighbor(const size_t id, const size_t currentSubs)
{
	for(size_t i = 0; i < atoms[currentSubs].bonds.size(); i++)
	{
		size_t neighborID = atoms[currentSubs].bonds[i].id;

		// if id is a nearest neighbor to one of current's partners,
		// id and current are secondNearestNeighbors
		if(isNearestNeighbor(id, atoms[neighborID-1]))
			return true;
	}
	return false;
}

/* Puts the atoms back in the periodic box
*/
void Model::keepInBox()
{        
	for(size_t i = 0; i < atoms.size(); i++)
	{
		while(atoms[i].coordinates.x < 0)
		{
				atoms[i].coordinates.x += BOX_SIZE[0];
				for(size_t j = 0; j < atoms[i].bonds.size(); j++)
					atoms[i].bonds[j].coordinates.x += BOX_SIZE[0];
		}
		while(atoms[i].coordinates.x >= BOX_SIZE[0])
		{
				atoms[i].coordinates.x -= BOX_SIZE[0];
				for(size_t j = 0; j < atoms[i].bonds.size(); j++)
					atoms[i].bonds[j].coordinates.x -= BOX_SIZE[0];
		}

		while(atoms[i].coordinates.y < 0)
		{
			atoms[i].coordinates.y += BOX_SIZE[1];
			for(size_t j = 0; j < atoms[i].bonds.size(); j++)
				atoms[i].bonds[j].coordinates.y += BOX_SIZE[1];
		}
		while(atoms[i].coordinates.y >= BOX_SIZE[1])
		{
			atoms[i].coordinates.y -= BOX_SIZE[1];
			for(size_t j = 0; j < atoms[i].bonds.size(); j++)
				atoms[i].bonds[j].coordinates.y -= BOX_SIZE[1];
		}

		while(atoms[i].coordinates.z < 0)
		{
			atoms[i].coordinates.z += BOX_SIZE[2];
			for(size_t j = 0; j < atoms[i].bonds.size(); j++)
				atoms[i].bonds[j].coordinates.z += BOX_SIZE[2];
		}
		while(atoms[i].coordinates.z >= BOX_SIZE[2])
		{
			atoms[i].coordinates.z -= BOX_SIZE[2];
			for(size_t j = 0; j < atoms[i].bonds.size(); j++)
				atoms[i].bonds[j].coordinates.z -= BOX_SIZE[2];
		}
	}
}

void Model::keepInInterfaceBox(double oldX)
{        
	for(size_t i = 0; i < atoms.size(); i++)
	{
		while(atoms[i].coordinates.x < 0)
		{
				atoms[i].coordinates.x += BOX_SIZE[0];
				for(size_t j = 0; j < atoms[i].bonds.size(); j++)
					atoms[i].bonds[j].coordinates.x += BOX_SIZE[0];
		}
		while(atoms[i].coordinates.x >= BOX_SIZE[0])
		{
				atoms[i].coordinates.x -= BOX_SIZE[0];
				for(size_t j = 0; j < atoms[i].bonds.size(); j++)
					atoms[i].bonds[j].coordinates.x -= BOX_SIZE[0];
		}

		while(atoms[i].coordinates.y < 0)
		{
			atoms[i].coordinates.y += BOX_SIZE[1];
			for(size_t j = 0; j < atoms[i].bonds.size(); j++)
				atoms[i].bonds[j].coordinates.y += BOX_SIZE[1];
		}
		while(atoms[i].coordinates.y >= BOX_SIZE[1])
		{
			atoms[i].coordinates.y -= BOX_SIZE[1];
			for(size_t j = 0; j < atoms[i].bonds.size(); j++)
				atoms[i].bonds[j].coordinates.y -= BOX_SIZE[1];
		}

		while(atoms[i].coordinates.z < 0)
		{
			atoms[i].coordinates.z += oldX;
			for(size_t j = 0; j < atoms[i].bonds.size(); j++)
				atoms[i].bonds[j].coordinates.z += oldX;
		}
		while(atoms[i].coordinates.z >= BOX_SIZE[2])
		{
			atoms[i].coordinates.z -= oldX;
			for(size_t j = 0; j < atoms[i].bonds.size(); j++)
				atoms[i].bonds[j].coordinates.z -= oldX;
		}
	}
}

/* Generates a crystal of SILICON or SIOX with sides of length BOX_SIZE*numBoxes and
8numBoxes^3 atoms
*/
void Model::makeCrystal(const int numBoxes, const int e, const bool interface)
{
	int id = 1;
	atoms.clear();

	BOX_SIZE[0] = 4 * B_0[e][e] / sqrt(3.0) * numBoxes; 
	//BOX_SIZE[0] = 4 * B_0[0][0] / sqrt(3.0) * numBoxes; 
	BOX_SIZE[2] = BOX_SIZE[1] = BOX_SIZE[0];

	for(int x = 0; x < numBoxes; x++)
	{
		for(int y = 0; y < numBoxes; y++)
		{
			for(int z = 0; z < numBoxes; z++)
			{
				atoms.push_back(Atom(Point(0+x,    0+y,    0+z   ), id++, e));
				atoms.push_back(Atom(Point(0.25+x, 0.25+y, 0.25+z), id++, e));
				atoms.push_back(Atom(Point(0.5+x,  0.5+y,  0+z   ), id++, e));
				atoms.push_back(Atom(Point(0.5+x,  0+y,    0.5+z ), id++, e));
				atoms.push_back(Atom(Point(0+x,    0.5+y,  0.5+z ), id++, e));
				atoms.push_back(Atom(Point(0.75+x, 0.75+y, 0.25+z), id++, e));
				atoms.push_back(Atom(Point(0.75+x, 0.25+y, 0.75+z), id++, e));
				atoms.push_back(Atom(Point(0.25+x, 0.75+y, 0.75+z), id++, e));
			}
		}
	}

	for(size_t i = 0; i < atoms.size(); i++)
	{
		atoms[i].coordinates.x *= (BOX_SIZE[0] / numBoxes);
		atoms[i].coordinates.y *= (BOX_SIZE[1] / numBoxes);
		atoms[i].coordinates.z *= (BOX_SIZE[2] / numBoxes);
	}
	
	if(!interface){
	    bondAtoms();
	}
}

void Model::insertH(const Point point, const int siId){

	int newId = atoms.size() + 1;
	Point p = closestPointWithPBC(point, atoms[siId-1].coordinates);
	atoms.push_back(Atom(point, newId, 6));
	double dist = p.distanceTo(atoms[newId-1].coordinates);
	
	atoms[newId-1].bonds.push_back(Partner(p, atoms[siId-1].id, dist));
	
	Vector displacement(p.x - atoms[newId-1].coordinates.x,
				p.y - atoms[newId-1].coordinates.y,
				p.z - atoms[newId-1].coordinates.z);

	atoms[siId-1].bonds.push_back(Partner(Point(atoms[siId-1].coordinates.x - displacement.x,
		        atoms[siId-1].coordinates.y - displacement.y,
				atoms[siId-1].coordinates.z - displacement.z),
				atoms[newId-1].id, dist));
	
	//atoms[siId-1].bonds.push_back(Partner(point, newId, dist));
	//atoms[newId-1].bonds.push_back(Partner(atoms[siId-1].coordinates, atoms[siId-1].id, dist));
	NUM_H++;

	keepInBox();
	
	return;
	
}
	
/* Makes a model with numAtoms silicon atoms with random positions
*/
// uses element 0
bool Model::makePreheatedCrystal(const int numAtoms)
{
	BOX_SIZE[0] = 2.715 * pow(numAtoms, 1/3.0);
	BOX_SIZE[2] = BOX_SIZE[1] = BOX_SIZE[0];
	while(true)
	{
		atoms.clear();
		for(int i = 0; i < numAtoms; i++)
			atoms.push_back(Atom(Point(getRandomDecimal()*BOX_SIZE[0], getRandomDecimal()*BOX_SIZE[1], getRandomDecimal()*BOX_SIZE[2]), i+1));

		if(bondAtoms())
			break;
		//cout << "failed\n";
	}

	initClosestNeighbors();
	for(int i = 0; i < 50 && relax() > 0; i++);
	keepInBox();

	return true;
}

bool Model::makeDefectedSystem()
{
	//Open output stream for error messages
	string errorFileName = FILE_NAME + ".error";
	ofstream fout;
	fout.open(errorFileName.c_str());
	
	BOX_SIZE[0] = 2.715 * pow(NEW_ATOMS, 1/3.0);
	BOX_SIZE[2] = BOX_SIZE[1] = BOX_SIZE[0];

	unsigned int numFail = 0;
	
	while(true)
	{
		if(numFail >= 100){
			fout << "Bonding failed after 100 attempts.  Program aborted." << endl;
			return false;
		}
		
		atoms.clear();
		int id = 1;
		
		for(int j = 0; j < N_TYPES; j++){

			for(int i = 0; i < nAtoms[j]; i++)
			{
				atoms.push_back(Atom(Point(getRandomDecimal()*BOX_SIZE[0], getRandomDecimal()*BOX_SIZE[1], getRandomDecimal()*BOX_SIZE[2]), id, nTypes[j]));
				id++;
			}
		}
		
		if(bondAtoms()){
			break;
		}
		else{
			numFail++;
		}
	}
	initClosestNeighbors();
	for(int i = 0; i < 50 && relax() > 0; i++);
	keepInBox();
	
	fout.close();

	return true;
}

bool Model::makePreheatedCrystal(const int numAtoms, const double x, const double y, const double z)
{
	BOX_SIZE[0] = x;
	BOX_SIZE[1] = y;
	BOX_SIZE[2] = z;
	
	while(true)
	{
		atoms.clear();
		for(int i = 0; i < numAtoms; i++)
			atoms.push_back(Atom(Point(getRandomDecimal()*BOX_SIZE[0], getRandomDecimal()*BOX_SIZE[1], getRandomDecimal()*BOX_SIZE[2]), i+1));

		if(bondAtoms())
			break;
	}

	initClosestNeighbors();
	for(int i = 0; i < 50 && relax() > 0; i++);
	keepInBox();

	return true;
}

void Model::deleteAtom(const int id)
{
	for(int i = 0; i < atoms[id-1].bonds.size(); i++)
	{
		int bondIndex = atoms[id-1].bonds[i].id-1;
		for (int j = 0; j < atoms[bondIndex].bonds.size(); j++)
		{
			if(atoms[bondIndex].bonds[j].id == id)
				atoms[bondIndex].bonds.erase(atoms[bondIndex].bonds.begin() + j);
		}
	}
	atoms.erase(atoms.begin() + id -1);

	for(int k = 0; k < atoms.size(); k++)
	{
		if(k > id)
			atoms[k].id--;
		for(int j =0; j < atoms[k].bonds.size(); j++)
		{
			if(atoms[k].bonds[j].id > id)
				atoms[k].bonds[j].id--;
		}
	}
}

void Model::deleteElement(const int type)
{
	for(int i = 0; i < atoms.size(); i++)
		if(atoms[i].type == type)
		{
			deleteAtom(i+1);
			i--;
		}
}

bool Model::makeInterface(string baseFile, const int type, const int numAtoms, bool fixed, bool ignore)
{
	/*if(!readInput(baseFile))
		return false;
	*/
	
	
	double oldBoxX = BOX_SIZE[2]; 
	double addition = 2.715 * pow(numAtoms, 1/3.0);
	BOX_SIZE[2] += addition;

	if(fixed)
		FIX_ATOMS = atoms.size();
	if(ignore)
		IGNORE_ATOMS = atoms.size();


	if(!randomPlacement(numAtoms, type, 0, 0, oldBoxX, BOX_SIZE[0], BOX_SIZE[1], addition))
		return false;
	
	//bondAtoms();

	return true;
}

bool Model::bondWithInterface(int bondTo)
{
	size_t n = atoms.size();

	// this avoids dead ends where atoms in less dense areas are not bonded
	//initAtomsInShell();
	//sort(atoms.begin(), atoms.end(), isMoreDesparate);

	// select an atom
	for(size_t i = 0; i < n; i++)
	{
		// while it doesn't have the needed bonds, bond
		while(atoms[i].bonds.size() < NUM_BONDS[atoms[i].type]) 
		{
			size_t closestSubscript = i;
			double closestDistance = INT_MAX;
			Point closestPoint = atoms[i].coordinates;

			// all elements 0-i already have max bonds and can be ignored
			for(size_t potentialPartnerSubs = i+1; 
				potentialPartnerSubs < n; 
				potentialPartnerSubs++)
			{
				if(atoms[potentialPartnerSubs].type == bondTo)
				{
				// a partner is valid if
				// 1. this partner does not have max bonds
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
	return true;
}

bool Model::randomPlacement(const int numAtoms, const int type, const double x, const double y, 
							const double z, const double x_length, const double y_length, const double z_length)
{
	OLD_X = z;
	int start = atoms.size();
	vector<Atom> temp = atoms;
	int j;

	for(j = 0; j < 100; j++)
	{ 
		atoms = temp;
		for(int i = 0; i < numAtoms; i++)
			atoms.push_back(Atom(Point((getRandomDecimal()*x_length) + x, (getRandomDecimal()*y_length) + y, (getRandomDecimal()*z_length) + z), start + i+1,type));

		if(bondWithInterface(type))
		{
			bool cont = true;
			float maxBond = BOX_SIZE[0] / 2;
			if(maxBond > BOX_SIZE[1]/2)
				maxBond = BOX_SIZE[1]/2;
			if(maxBond > BOX_SIZE[2]/2)
				maxBond = BOX_SIZE[2]/2;

			for(int i = 0; i < atoms.size(); i++)
			{
				for(int j = 0; j < atoms[i].bonds.size(); j++)
				{
					if(atoms[i].bonds[j].dist >= maxBond)
						cont = false;
				}
			}

			if(cont)
				break;
		}
	}
	
	if(j >= 100)
		return false;
	else
	{
		initClosestNeighbors();
		for(int i = 0; i < 50 && relax() > 0; i++);
		if(INTER==true){
			keepInInterfaceBox(OLD_X);
		}
		else{
			keepInBox();
		}
		return true;
	}
    
}

/* Performs numSwaps annealing transpositions
Sends energy information to a file with a given filename, if one is specified
Otherwise, there is no output, which is much faster
NOTE: Set kT in input file before calling this
*/
void Model::performAnnealing(const int numSwaps, const int numVSwaps, string filename)
{
	double kT0 = kT;     // initial temperature
	double kTf = 0.1;     // final temperature
	bool OUTPUT_ON = true;
	ofstream fout;
	int acceptedSwaps = 0;

	if(OUTPUT_ON)
	{
		fout.open(filename.c_str());
	         fout << "Swap #, initial energy, "
			<< "# step, final energy, "
			<< "random number, kT, "
			<< "Metropolis number" << endl;
}

	double initEnergy = getEnergy();
	//The annealing function needs to save the lowest energy configuration, which is not necessarily the final one
	double minEnergy = initEnergy;
	
	vector<Atom> minAtoms = atoms;
	
	for(int i = 0; i < numSwaps; i++)
	{
		//cout << "Swap: " << i << endl;
		vector<Atom> backup = atoms;
		double initEnergy = getEnergy();
		randomSwap();
		int numSteps = relax();
		if(INTER==true){
			keepInInterfaceBox(OLD_X);
		}
		else{
			keepInBox();
		}

		double finEnergy = getEnergy();
		double kT_local = kT*finEnergy* pow(1- (double)i/(double)numSwaps, 4)/minEnergy;
		double rand = getRandomDecimal();
                double metProb = getMetropolisProb(initEnergy, finEnergy, kT_local);
		string font = "<font color=blue>";

		if(rand >= metProb)
		{
    		// rejected → do nothing except restore atoms
    			atoms = backup;
		}		
		else{
    		// accepted
    			acceptedSwaps++;
 			//cout << "Swap: " << i << endl;
    			if(RELAX_V == 1 && i >= numVSwaps){
        			relaxVolume();
        			finEnergy = getEnergy();
    			}

    			if(finEnergy < minEnergy){
        			minEnergy = finEnergy;
        			minAtoms = atoms;
    			}		

    		// ✅ WRITE ONLY ACCEPTED SWAPS
    			if(OUTPUT_ON){ 
        			fout <<  i+1 << ", "
             			<< initEnergy << ", " 
             			<< numSteps << ", "
             			<< finEnergy << ", "
             			<< rand << ", "
             		        << kT_local << ", "	
				<< metProb << endl;
    			}	
		}
//	atoms = minAtoms;
	}
	if(RELAX_V == 1){
		relaxVolume();
		minEnergy = getEnergy();
	}
	if(OUTPUT_ON){ 
    		cout << "Total Number of accepted swaps: " 
         	<< acceptedSwaps << "";

    		cout << "Minimum energy: " 
         	<< minEnergy << "";
		fout.close();
	} 
}

/* Performs a random WWW bond transposition on the system
*/
void Model::randomSwap()
{
        size_t  aSubs, bSubs, cSubs, dSubs, initR, trialR;
        unsigned int randNum;
        size_t randIndex;
        unsigned int maxPercent = atomPercent[0];
        unsigned int swapType;



	while(true)
	{
		// choose A based on atom swap probabilities
		randNum = getRandomInteger();
		if(N_TYPES > 1){
			
			for(int i = 0; i < atomPercent.size(); i++){
				
				nIds.clear(); 
				
				if(randNum <= maxPercent){
					
					swapType = nTypes[i];
					
					for(int j = 0; j < atoms.size(); j++){
						
						if(atoms[j].type == swapType){
							
							nIds.push_back(atoms[j].id);
							
						}
					}
					
					randIndex = (rand() % (nIds.size()));
					
					aSubs = nIds[randIndex] - 1;
					

					break;
					
				}
				else if(i < atomPercent.size()-1){
					
					maxPercent = maxPercent + atomPercent[i+1];
				}
				/*else{
					aSubs = (rand() % (atoms.size()));
				}*/
			}
		}
		else{
			aSubs = (rand() % (atoms.size()));
			
		}

		// choose B
		initR = rand() % atoms[aSubs].bonds.size();
		bSubs = atoms[aSubs].bonds[initR].id - 1;

		if(bSubs < IGNORE_ATOMS)
			continue;

		// choose C
		initR = rand() % atoms[bSubs].bonds.size();
		cSubs = atoms[bSubs].bonds[initR].id - 1;

		if(cSubs < IGNORE_ATOMS)
			continue;

		// there is a problem with C
		if(cSubs == aSubs || isNearestNeighbor(cSubs+1, atoms[aSubs]))
		{
			// so try another C
			trialR = (initR + 1) % atoms[bSubs].bonds.size();
			cSubs = atoms[bSubs].bonds[trialR].id - 1;

			while((cSubs == aSubs || isNearestNeighbor(cSubs+1, atoms[aSubs])) && trialR != initR)
			{
				trialR = (trialR + 1) % atoms[bSubs].bonds.size();
				cSubs = atoms[bSubs].bonds[trialR].id - 1;
			}

			// no value of C will work
			if(trialR == initR)
				continue;
			if(cSubs < IGNORE_ATOMS)
				continue;
		}

		// choose D
		initR = rand() % atoms[cSubs].bonds.size();
		dSubs = atoms[cSubs].bonds[initR].id - 1;

		if(dSubs < IGNORE_ATOMS)
			continue;

		// there is a problem with D
		if(dSubs == bSubs || isNearestNeighbor(dSubs+1, atoms[bSubs]))
		{
			// so try another D
			trialR = (initR + 1) % atoms[cSubs].bonds.size();
			dSubs = atoms[cSubs].bonds[trialR].id - 1;

			while((dSubs == bSubs || isNearestNeighbor(dSubs+1, atoms[bSubs])) && trialR != initR)
			{
				trialR = (trialR + 1) % atoms[cSubs].bonds.size();
				dSubs = atoms[cSubs].bonds[trialR].id - 1;
			}

			// no value of D will work
			if(trialR == initR)
				continue;
			if(dSubs < IGNORE_ATOMS)
				continue;
		}

		break;
	}

	// swapping them to be a-c-b-d
	//cout << "Swapping " << aSubs+1 << " " << bSubs+1 << " " << cSubs+1 << " " << dSubs+1 << endl;
	swapBonds(aSubs, bSubs, cSubs, dSubs);
}

/*
Reads Constants From a Model File
*/

void Model::readConstants(const char* filename) 
{	
	char ws;
	//Open file for reading constants
	ifstream fin(filename); 
	fin >> NUM_ELEMENTS;
	fin.get(ws);
	while(!fin.eof())
	{
		char c;
		string constName;

		while(fin >> c)
		{   
			if(c == '*')
			{
				fin >> constName;

				if(constName == "NAME")
				{
					string defaultV, elName;
					unsigned int element;
					fin >> defaultV;
					NAME = vector<string>(NUM_ELEMENTS, defaultV);
					fin.get(ws);
					while(!fin.eof() && fin.peek() != '*')
					{

						fin >> element >> elName;
						NAME[element] = elName;
						while(isspace(fin.peek())){
							fin.get(ws);
						}
					}
				}
			        else if(constName == "NUM_BONDS")
				{
					unsigned int defaultV, element, numBonds;
					fin >> defaultV;
					NUM_BONDS = vector<unsigned int>(NUM_ELEMENTS, defaultV);
					fin.get(ws);
					while(!fin.eof() && fin.peek() !=  '*') 
					{
						fin >> element >> numBonds;
						NUM_BONDS[element] = numBonds;
						while(isspace(fin.peek())){
							fin.get(ws);
						}
					}
				}
				else if(constName == "K_OMEGA")
				{	
					double defaultV, k_Omega;
					fin >> defaultV;
					unsigned int element1, element2, element3;

					K_OMEGA = vector<vector<vector<double> > >(NUM_ELEMENTS, vector<vector<double> >(NUM_ELEMENTS, vector<double>(NUM_ELEMENTS, defaultV)));
					fin.get(ws);
					while(!fin.eof() && fin.peek() !=  '*') 
					{
						fin >> element1 >> element2 >> element3 >> k_Omega;
						K_OMEGA[element1][element2][element3] = k_Omega;
						K_OMEGA[element1][element3][element2] = k_Omega;
						while(isspace(fin.peek())){
							fin.get(ws);
						}
					}
				}
				else if(constName == "COS_ANGLE")
				{
					double defaultV, cosAngle;
					unsigned int element;
					fin >> defaultV;
					COS_ANGLE = vector<double>(NUM_ELEMENTS, defaultV);
					fin.get(ws);
					while(!fin.eof() && fin.peek() !=  '*') 
					{
						fin >> element >> cosAngle;
						COS_ANGLE[element] = cosAngle;
						while(isspace(fin.peek())){
							fin.get(ws);
						}
					}
				}
				else if(constName == "K_B")
				{
					double defaultV, k_B;
					unsigned int element1, element2;
					fin >> defaultV;
					K_B = vector<vector<double> >(NUM_ELEMENTS, vector<double>(NUM_ELEMENTS, defaultV));
					fin.get(ws);
					while(!fin.eof() && fin.peek() !=  '*') 
					{
						fin >> element1 >> element2 >> k_B;
						K_B[element1][element2] = k_B;
						K_B[element2][element1] = k_B;
						while(isspace(fin.peek())){
							fin.get(ws);
						}
					}
				}
				else if(constName == "B_0")
				{
					double defaultV, b_0;
					unsigned int element1, element2;
					fin >> defaultV;
					B_0 = vector<vector<double> >(NUM_ELEMENTS, vector<double>(NUM_ELEMENTS, defaultV));
					fin.get(ws);
					while(!fin.eof() && fin.peek() !=  '*') 
					{
						fin >> element1 >> element2 >> b_0;
						B_0[element1][element2] = b_0;
						B_0[element2][element1] = b_0;
						while(isspace(fin.peek())){
							fin.get(ws);
						}	
					}
				}
				else if(constName == "D_0")
				{ 
					double defaultV;
					double d_0;
					unsigned int element1, element2;

					fin >> defaultV;
					D_0 = vector<vector<double> >(NUM_ELEMENTS, vector<double>(NUM_ELEMENTS, defaultV));
					fin.get(ws);
					while(!fin.eof() && fin.peek() !=  '*') 
					{
						fin >> element1 >> element2 >> d_0;
						D_0[element1][element2] = d_0;
						D_0[element2][element1] = d_0;
						while(isspace(fin.peek())){
							fin.get(ws);
						}
					}
				}
				else if(constName == "GAMMA")
				{ 
					double defaultV;
					double gamma;
					unsigned int element1, element2;

					fin >> defaultV;
					GAMMA = vector<vector<double> >(NUM_ELEMENTS, vector<double>(NUM_ELEMENTS, defaultV));
					fin.get(ws);
					while(!fin.eof() && fin.peek() !=  '*') 
					{
						fin >> element1 >> element2 >> gamma;
						GAMMA[element1][element2] = gamma;
						GAMMA[element2][element1] = gamma;
						while(isspace(fin.peek())){
							fin.get(ws);
						}
					}
				}
				else if(constName == "MASS")
				{
					double defaultV;
					double mass;
					unsigned int element;
					fin >> defaultV;
					MASS = vector<double>(NUM_ELEMENTS, defaultV);
					fin.get(ws);
					while(!fin.eof() && fin.peek() !=  '*') 
					{
						fin >> element >> mass;
						MASS[element] = mass;
						while(isspace(fin.peek())){
							fin.get(ws);
						}
					}
				}
			}
		}
	}
	
}

bool Model::readInput()
{
	char ws;
	string strin;
	unsigned int in = 0;
	atoms.clear();
	
	ifstream fin;
	fin.open(C_FILE_NAME.c_str());
	
	bool error = false;
	
	fin >> FILE_NAME;
	fin >> P_BONDS >> P_VASP >> WRITE_NN;
	while(isspace(fin.peek())){
		fin.get(ws);
	}
	getline(fin, strin);
	istringstream iss(strin);
	while(iss >> in){
		nTypes.push_back(in);
	}
	N_TYPES = nTypes.size();

	getline(fin, strin);
	istringstream niss(strin);
	while(niss >> in){
		atomPercent.push_back(in);
	}
	
	getline(fin, strin);
	istringstream siss(strin);
	while(siss >> in){
		nAtoms.push_back(in);
		NEW_ATOMS = NEW_ATOMS + in;
	}
	
	fin >> kT;
	fin >> NUM_SWITCHES;
	fin >> BOX_SIZE[0] >> BOX_SIZE[1] >> BOX_SIZE[2];
	fin >> RELAX_V >> RELAX_V_SWITCHES;

	//Open output stream for error messages
	string errorFileName = FILE_NAME + ".error";
	ofstream fout;
	fout.open(errorFileName.c_str());
	
	if(P_BONDS != 0 && P_BONDS != 1){
		fout << "The leftmost number on line 2 (.b file creation) in the file " << C_FILE_NAME << " should be 0 or 1.  Program aborted." << endl;
		error = true;
	}
	
	if(P_VASP != 0 && P_VASP != 1){
		fout << "The middle number on line 2 (.vasp file creation) in the file " << C_FILE_NAME << " should be 0 or 1.  Program aborted." << endl;
		error = true;
	}

	if(WRITE_NN != 0 && WRITE_NN != 1){
		fout << "The rightmost number on line 2 (.nn file creation) in the file " << C_FILE_NAME << " should be 0 or 1.  Program aborted." << endl;
		error = true;
	}

	for(int j = 0; j < N_TYPES; j++){
		if(nTypes[j] >= NUM_ELEMENTS){
			fout << "You have entered an invalid atom type on line 3.  Program aborted." << endl;
			error = true;	
		}
		
	}

	if(N_TYPES == 1 && nTypes[0] == 2){
		fout << "You have entered only H atoms on line 3.  You cannot have only H atoms in a model.  Program aborted." << endl;
		error = true;
	}

	for(int j = 0; j < N_TYPES; j++){
		if(nAtoms[j] == 0){
			fout << "You entered 0 atoms for one or more atom types on line 5.  Please remove the corresponding atom type(s) from line 3.  Program aborted." << endl;
			error = true;
		}
	}
	
	unsigned int totalPercent = 0;
	for(int j = 0; j < atomPercent.size(); j++){
		totalPercent = totalPercent + atomPercent[j];
	}

	if(totalPercent != 100){
		fout << "The numbers on line 4 should add up to 100.  Program aborted." << endl;
		error = true;
	}

	if(N_TYPES != nAtoms.size()){
		fout << "The number of atom types on line 3 does not match the numbers of atoms on line 5.  Program aborted." << endl;
		error = true;
	}
	
	if(kT < 0 || kT > 20){
		fout << "The kT value (line 6) in the file " << C_FILE_NAME << " needs to be between 0 and 20.  Program aborted." << endl;
		error = true;
	}
	
	if (static_cast<int>(NUM_SWITCHES) == NUM_SWITCHES){
		I_NUM_SWITCHES = static_cast<int>(NUM_SWITCHES);
	}
	else{
		fout << "The number of bond swaps (line 7) in the file " << C_FILE_NAME << " needs to be an integer value.  Program aborted." << endl;
		error = true;
	}
	
	if(RELAX_V != 0 && RELAX_V != 1){
		fout << "The leftmost number on line 9 in the file " << C_FILE_NAME << " should be 0 or 1.  Program aborted." << endl;
		error = true;
	}
	
	if (static_cast<int>(RELAX_V_SWITCHES) == RELAX_V_SWITCHES){
		I_RELAX_V_SWITCHES = static_cast<int>(RELAX_V_SWITCHES);
	}
	else{
		fout << "The number of bond swaps before volume relaxation (line 9) in the file " << C_FILE_NAME << " needs to be an integer value.  Program aborted." << endl;
		error = true;
	}
	
	if(I_RELAX_V_SWITCHES > I_NUM_SWITCHES){
		fout << "The number of bond swaps before volume relaxation (line 9) in the file " << C_FILE_NAME << " cannot be greater than the total number of bond swaps performed.  Program aborted." << endl;
		error = true;
	}
	
	//if(error == true)
	//	return false;
		
	int id, numBonds, type; 
	double x, y, z;

	if(error == false){
		while(fin >> id >> type >> x >> y >> z >> numBonds){
			/*if(x > BOX_SIZE[0])
				x = x - BOX_SIZE[0];
			else if(x < 0)
				x = BOX_SIZE[0] + x;
			if(y > BOX_SIZE[1])
				y = y - BOX_SIZE[1];
			else if(y < 0)
				y = BOX_SIZE[1] + y;
			if(z > BOX_SIZE[2])
				z = z - BOX_SIZE[2];
			else if(z < 0)
				z = BOX_SIZE[2] + z;
			*/

			/* ClosestPointWithPBC(const Point &p1, const Point &p2) */

			atoms.push_back(Atom(Point(x, y, z), id));
			atoms[atoms.size()-1].type = type;
			
			Point temp=Point(x,y,z);	//test					

			for(int i = 0; i < numBonds; i++) 
			{
				fin >> id >> x >> y >> z;
				Point compTemp=Point(x,y,z); //test
				atoms[atoms.size()-1].bonds.push_back(Partner(closestPointWithPBC(temp, compTemp), id, atoms[atoms.size()-1].coordinates.distanceTo(closestPointWithPBC(temp, compTemp))));
			}
		}
	}

	fin.close();
	
	if(error == true)
		return false;
		
	fout.close();

	return true;
}

/* Reads a model from a file
*/
void Model::readModel()
{
	atoms.clear(); 
	//Remove existing atoms 
	//ifstream fin(filename.c_str()); 

	//Store length of box 
	cin >> BOX_SIZE[0] >> BOX_SIZE[1] >> BOX_SIZE[2];
	
	int id, numBonds, type; 
	double x, y, z; 

	while(cin >> id >> type >> x >> y >> z >> numBonds) //While line is not empty 
	{
		atoms.push_back(Atom(Point(x, y, z), id));
		atoms[atoms.size()-1].type = type;

		for(int i = 0; i < numBonds; i++) 
		{
			cin >> id >> x >> y >> z;
			atoms[atoms.size()-1].bonds.push_back(Partner(Point(x, y, z), id, atoms[atoms.size()-1].coordinates.distanceTo(Point(x,y,z))));
		}
	}
}

/* Allows the atoms to move until the forces become very small
*/
int Model::relax()
{
	int numSteps;
	for(numSteps = 0; ; numSteps++)
	{
		setForces();
		setAccelerations();
		double biggestForce = getLargestForce();

		if(biggestForce < 0.1 || numSteps >= 300)
			break;
		else if(biggestForce < 0.5)
			setPositions(15);
		else if(biggestForce < 2)
			setPositions(10);	
		else if(biggestForce < 10)
			setPositions(7);
		else if(biggestForce < 20)
			setPositions(5);
		else
			setPositions(2);
	
		
	/*	if(INTER==true){
			keepInInterfaceBox(OLD_X);
		}
		else{
			keepInBox();
		}
	*/
	}

	return numSteps;
}

void Model::relaxVolume(){

	double percentChange = 0.05;

	double firstSize0 = BOX_SIZE[0];
	double firstSize1 = BOX_SIZE[1];
	double firstSize2 = BOX_SIZE[2];
	
	double origEnergy = getEnergy();
	double newEnergy = 0.0;
	
	changeBoxSize((BOX_SIZE[0])*(1.00+percentChange));
	
	newEnergy = getEnergy();
	
	if(newEnergy < origEnergy){

		changeBoxSize(firstSize0);
	
		while(newEnergy < origEnergy){
		
			changeBoxSize((BOX_SIZE[0])*(1.00+percentChange));
			
			newEnergy = getEnergy();
				
			if(newEnergy < origEnergy){
				origEnergy = newEnergy;
				firstSize0 = BOX_SIZE[0];
				firstSize1 = BOX_SIZE[1];
				firstSize2 = BOX_SIZE[2];
			}
			else{
				changeBoxSize(firstSize0);
				break;
			}
		}
	}
	
	else if(newEnergy > origEnergy){
	
		changeBoxSize((BOX_SIZE[0])*(1.00-percentChange));
		
		newEnergy = getEnergy();
		
		if(newEnergy < origEnergy){
		
			firstSize0 = BOX_SIZE[0];
			firstSize1 = BOX_SIZE[1];
			firstSize2 = BOX_SIZE[2];
		
			while(newEnergy < origEnergy){
			
				changeBoxSize((BOX_SIZE[0])*(1.00-percentChange));
				
				newEnergy = getEnergy();
					
				if(newEnergy < origEnergy){
					origEnergy = newEnergy;
					firstSize0 = BOX_SIZE[0];
					firstSize1 = BOX_SIZE[1];
					firstSize2 = BOX_SIZE[2];
				}
				else{
					changeBoxSize(firstSize0);
					break;
				}
			}
		}
    }
	else{
		changeBoxSize(BOX_SIZE[0]);
    }

	return;
}

/* Sets the accelerations based on the forces
*/
void Model::setAccelerations()
{
	size_t n = atoms.size();

	for(size_t i = FIX_ATOMS; i < n; i++)
	{
		// converts the force to kg A/fs^2, then divides by the atom mass in kg
		atoms[i].acceleration.x = atoms[i].force.x * FORCE_CONVERSION / MASS[atoms[i].type];
		atoms[i].acceleration.y = atoms[i].force.y * FORCE_CONVERSION / MASS[atoms[i].type];
		atoms[i].acceleration.z = atoms[i].force.z * FORCE_CONVERSION / MASS[atoms[i].type];
	}
}

/* Sets the angular energy of the system
*/
void Model::setAngularEnergy()
{
	/* angular energy belongs to the center atom */
	static size_t i,j,k,n=atoms.size(),jSubs,kSubs;
	static Point pj,pk;

	// select atom
	for( i = 0; i < n; i++)
	{
		//if(SHELL_ONLY && atoms[i].shellNumber > CLUSTER_MAX_SHELL+1)	continue;

		atoms[i].angularEnergy = 0;
		// select first bond
		for( j = 0; j < atoms[i].bonds.size(); j++)
		{
			jSubs=atoms[i].bonds[j].id -1;
			pj = closestPointWithPBC(atoms[i].coordinates,atoms[jSubs].coordinates);

			// select second bond
			for(size_t k = j+1; k < atoms[i].bonds.size(); k++)
			{
				kSubs = atoms[i].bonds[k].id-1;
				pk=closestPointWithPBC(atoms[i].coordinates, atoms[kSubs].coordinates);
				// 0.5 0.647 2.35^2 (cosijk + 1/3)^2
				atoms[i].angularEnergy += (K_OMEGA[atoms[i].type][atoms[jSubs].type][atoms[kSubs].type] * 
					pow(cosijk(atoms[i].coordinates, pj,pk,atoms[i].bonds[j].dist,atoms[i].bonds[k].dist) - COS_ANGLE[atoms[i].type], 2) / 2);
			}
		}
	}
}

/* Sets the total energy of the system
*/
void Model::setEnergy()
{
	if(RADIAL_ON)	setRadialEnergy();
	if(ANGULAR_ON)	setAngularEnergy();
	if(REPEL_ON)	setRepulsiveEnergy();
}

/* Sets the forces on the atoms
*/
void Model::setForces()
{
	size_t n = atoms.size();
	Partner temp;

	size_t i, j, k, jSubs, kSubs;
	double multiplier, tempNum, cosijk, tempJ, tempK, dist;

	Point p;

	// reset all forces
	for(i = 0; i < n; i++)
	{
		if((SHELL_ONLY && atoms[i].shellNumber > CLUSTER_MAX_SHELL))	
			//this line was changed as a test
			continue;

		atoms[i].force.x = 0;
		atoms[i].force.y = 0;
		atoms[i].force.z = 0;
	}

	for(i = 0; i < n; i++)
	{
		if((SHELL_ONLY && atoms[i].shellNumber > CLUSTER_MAX_SHELL))	
			continue;

		for(j = 0; j < atoms[i].bonds.size(); j++)
		{
			jSubs = atoms[i].bonds[j].id-1;

			// set radial forces of two partner atoms at a time
			// twice as efficient, since it takes the "equal and opposite" rule into account
			if(RADIAL_ON  &&  i < jSubs)
			{
				multiplier = (K_B[atoms[i].type][atoms[atoms[i].bonds[j].id-1].type] * (atoms[i].bonds[j].dist-B_0[atoms[i].type][atoms[atoms[i].bonds[j].id-1].type]) / atoms[i].bonds[j].dist);

				tempNum = multiplier * (atoms[i].bonds[j].coordinates.x - atoms[i].coordinates.x);
				atoms[i].force.x += tempNum;
				atoms[jSubs].force.x -= tempNum;

				tempNum = multiplier * (atoms[i].bonds[j].coordinates.y - atoms[i].coordinates.y);
				atoms[i].force.y += tempNum;
				atoms[jSubs].force.y -= tempNum;

				tempNum = multiplier * (atoms[i].bonds[j].coordinates.z - atoms[i].coordinates.z);
				atoms[i].force.z += tempNum;
				atoms[jSubs].force.z -= tempNum;
			} // end if

			// set angular force of three atoms at a time
			// i is the "middle" atom and j/k are endpoints of the angle
			for(k = j+1; k < atoms[i].bonds.size() && ANGULAR_ON; k++)
			{
				kSubs = atoms[i].bonds[k].id-1;

				// avoids recomputing these repeatedly
				cosijk = ((atoms[i].bonds[j].coordinates.x - atoms[i].coordinates.x) * (atoms[i].bonds[k].coordinates.x - atoms[i].coordinates.x) + 
					(atoms[i].bonds[j].coordinates.y - atoms[i].coordinates.y) * (atoms[i].bonds[k].coordinates.y - atoms[i].coordinates.y) + 
					(atoms[i].bonds[j].coordinates.z - atoms[i].coordinates.z) * (atoms[i].bonds[k].coordinates.z - atoms[i].coordinates.z)) / 
					(atoms[i].bonds[j].dist * atoms[i].bonds[k].dist);
				multiplier = K_OMEGA[atoms[i].type][atoms[atoms[i].bonds[j].id-1].type][atoms[atoms[i].bonds[k].id-1].type] * 
					(COS_ANGLE[atoms[i].type]-cosijk);

				// find force on j
				tempJ = multiplier * ((atoms[i].bonds[k].coordinates.x-atoms[i].coordinates.x) / atoms[i].bonds[k].dist +
					(atoms[i].coordinates.x-atoms[i].bonds[j].coordinates.x) * cosijk / atoms[i].bonds[j].dist) / atoms[i].bonds[j].dist;
				atoms[jSubs].force.x += tempJ;
				atoms[i].force.x -= tempJ;

				tempJ = multiplier * ((atoms[i].bonds[k].coordinates.y-atoms[i].coordinates.y) / atoms[i].bonds[k].dist +
					(atoms[i].coordinates.y-atoms[i].bonds[j].coordinates.y) * cosijk / atoms[i].bonds[j].dist) / atoms[i].bonds[j].dist;
				atoms[jSubs].force.y += tempJ;
				atoms[i].force.y -= tempJ;

				tempJ = multiplier * ((atoms[i].bonds[k].coordinates.z-atoms[i].coordinates.z) / atoms[i].bonds[k].dist +
					(atoms[i].coordinates.z-atoms[i].bonds[j].coordinates.z) * cosijk / atoms[i].bonds[j].dist) / atoms[i].bonds[j].dist;
				atoms[jSubs].force.z += tempJ;
				atoms[i].force.z -= tempJ;

				// find force on k
				tempK = multiplier * ((atoms[i].bonds[j].coordinates.x-atoms[i].coordinates.x) / atoms[i].bonds[j].dist +
					(atoms[i].coordinates.x-atoms[i].bonds[k].coordinates.x) * cosijk / atoms[i].bonds[k].dist) / atoms[i].bonds[k].dist;
				atoms[kSubs].force.x += tempK;
				atoms[i].force.x -= tempK;

				tempK = multiplier * ((atoms[i].bonds[j].coordinates.y-atoms[i].coordinates.y) / atoms[i].bonds[j].dist +
					(atoms[i].coordinates.y-atoms[i].bonds[k].coordinates.y) * cosijk / atoms[i].bonds[k].dist) / atoms[i].bonds[k].dist;
				atoms[kSubs].force.y += tempK;
				atoms[i].force.y -= tempK;

				tempK = multiplier * ((atoms[i].bonds[j].coordinates.z-atoms[i].coordinates.z) / atoms[i].bonds[j].dist +
					(atoms[i].coordinates.z-atoms[i].bonds[k].coordinates.z) * cosijk / atoms[i].bonds[k].dist) / atoms[i].bonds[k].dist;
				atoms[kSubs].force.z += tempK;
				atoms[i].force.z -= tempK;
			} // end k loop
		} // end j loop

		if(REPEL_ON)
		{
			// repulsive force term
			for(j = 0; j < atoms[i].closestNeighbors.size(); j++)
			{
				// find closest version of atom
				p = closestPointWithPBC(atoms[i].coordinates, atoms[atoms[i].closestNeighbors[j]-1].coordinates);
				dist = p.distanceTo(atoms[i].coordinates);

				// if the distance is greater than D_0, there is no force
				if(dist < D_0[atoms[i].type][atoms[atoms[i].closestNeighbors[j]-1].type])
				{
					multiplier = 3 * GAMMA[atoms[i].type][atoms[atoms[i].closestNeighbors[j]-1].type] * pow(D_0[atoms[i].type][atoms[atoms[i].closestNeighbors[j]-1].type] - dist, 2) / dist;
					atoms[i].force.x += (multiplier * (atoms[i].coordinates.x - p.x));
					atoms[i].force.y += (multiplier * (atoms[i].coordinates.y - p.y));
					atoms[i].force.z += (multiplier * (atoms[i].coordinates.z - p.z));
				}
			}
		} // end if repel on 
	} // end i loop

}

/* Sets the positions of the atoms based on the acceleration
*/
void Model::setPositions(const double timestep)
{
	Vector displacement;
	double multiplier = pow(timestep, 2) / 2.0;
	size_t jSubs;

	size_t n = atoms.size();

	for(size_t i = FIX_ATOMS; i < n; i++)
	{
		if((SHELL_ONLY && atoms[i].shellNumber > CLUSTER_MAX_SHELL) )
			//|| (i < FIX_ATOMS))	
			continue;

		// calculate displacement
		displacement.x = atoms[i].acceleration.x * multiplier;
		displacement.y = atoms[i].acceleration.y * multiplier;
		displacement.z = atoms[i].acceleration.z * multiplier;
		// adjust coordinates
		atoms[i].coordinates.x += displacement.x;
		atoms[i].coordinates.y += displacement.y;
		atoms[i].coordinates.z += displacement.z;

		// move the bonds
		for(size_t j = 0; j < atoms[i].bonds.size(); j++)
		{
			// adjust bond lengths
			atoms[i].bonds[j].dist = atoms[i].coordinates.distanceTo(atoms[i].bonds[j].coordinates);

			// moves partner atoms
			jSubs = atoms[i].bonds[j].id - 1;
			for(size_t k = 0; k < atoms[jSubs].bonds.size(); k++)
			{
				if(atoms[i].id == atoms[jSubs].bonds[k].id)
				{
					atoms[jSubs].bonds[k].coordinates.x += displacement.x;
					atoms[jSubs].bonds[k].coordinates.y += displacement.y;
					atoms[jSubs].bonds[k].coordinates.z += displacement.z;
					atoms[jSubs].bonds[k].dist = atoms[i].bonds[j].dist;
				}
			}
		}
	}
}

/* Sets the radial energy of the atoms
*/
void Model::setRadialEnergy()
{
	/* radial energy is divided between the two atoms */
	static size_t i,j,n=atoms.size(),jSubs;

	// select atom
	for( i = 0; i < n; i++)
	{

		//if(SHELL_ONLY && atoms[i].shellNumber > CLUSTER_MAX_SHELL+1)	continue;

		atoms[i].radialEnergy = 0;
		// select first bond
		for( j = 0; j < atoms[i].bonds.size(); j++)
		{
//   BRT - correction to jSubs 06/2014
			jSubs = atoms[i].bonds[j].id - 1;
			// divided by 4 not 2 because two atoms share the energy
			atoms[i].radialEnergy += K_B[atoms[i].type][atoms[jSubs].type] * pow(atoms[i].bonds[j].dist - B_0[atoms[i].type][atoms[jSubs].type], 2) / 4;
		}
	}
}

/* Sets the repulsive energy of the atoms
*/
void Model::setRepulsiveEnergy()
{
	/* repulsive energy is divided between two atoms */

	static size_t i,j,n=atoms.size();
	static Point tempPt;
	static double tempDist;	

	for( i = 0; i < n; i++)	// subscript 1
	{
		//if(SHELL_ONLY && atoms[i].shellNumber > CLUSTER_MAX_SHELL+1)	continue;

		atoms[i].repulsiveEnergy = 0;
		for( j = 0; j < atoms[i].closestNeighbors.size(); j++)	
		{
			// find closest version of this point that might be repelling
			tempPt = closestPointWithPBC(atoms[i].coordinates, atoms[atoms[i].closestNeighbors[j]].coordinates);
			tempDist = tempPt.distanceTo(atoms[i].coordinates);

			if(tempDist < D_0[atoms[i].type][atoms[atoms[i].closestNeighbors[j]].type])				// add term to energy ~ 0.5 (3.84 - |rm - rn|)^3
				atoms[i].repulsiveEnergy += (GAMMA[atoms[i].type][atoms[atoms[i].closestNeighbors[j]].type] * pow(D_0[atoms[i].type][atoms[atoms[i].closestNeighbors[j]].type] - tempDist, 3) / 2);
		}
	}
}

/* Swaps the bonds of four atoms based on the WWW algorithm
*/
void Model::swapBonds(const size_t aSubs, const size_t bSubs, 
					  const size_t cSubs, const size_t dSubs)
{

	for(size_t i = 0; i < 4; i++)
	{
		Point temp;
		// a-b becomes a-c
		if(i < atoms[aSubs].bonds.size() && atoms[aSubs].bonds[i].id == bSubs+1)
		{
			temp = closestPointWithPBC(atoms[aSubs].coordinates, atoms[cSubs].coordinates);
			atoms[aSubs].bonds[i] = Partner(temp, cSubs+1, temp.distanceTo(atoms[aSubs].coordinates)); 
		}

		// b-a becomes b-d
		if(i < atoms[bSubs].bonds.size() && atoms[bSubs].bonds[i].id == aSubs+1)
		{
			temp = closestPointWithPBC(atoms[bSubs].coordinates, atoms[dSubs].coordinates);
			atoms[bSubs].bonds[i] = Partner(temp, dSubs+1, temp.distanceTo(atoms[bSubs].coordinates));
		}

		// c-d becomes c-a
		if(i < atoms[cSubs].bonds.size() && atoms[cSubs].bonds[i].id == dSubs+1)
		{
			temp = closestPointWithPBC(atoms[cSubs].coordinates, atoms[aSubs].coordinates);
			atoms[cSubs].bonds[i] = Partner(temp, aSubs+1, temp.distanceTo(atoms[cSubs].coordinates));
		}

		// d-c becomes d-b
		if(i < atoms[dSubs].bonds.size() && atoms[dSubs].bonds[i].id == cSubs+1)
		{
			temp = closestPointWithPBC(atoms[dSubs].coordinates, atoms[bSubs].coordinates);
			atoms[dSubs].bonds[i] = Partner(temp, bSubs+1, temp.distanceTo(atoms[dSubs].coordinates));
		}

	}

	initClosestNeighbors();
	if(SHELL_ONLY)	initShellNumbers(aSubs, bSubs, cSubs, dSubs);
}

/* Sends the bond information to a file
*/
void Model::writeBonds(string filename)
{
	ofstream fout(filename.c_str());

	fout << "<html><head><title>Atom Bonds</title></head>\n"
		<< "<body><table border=1>\n"
		<< "<tr><td>Atom ID</td><td>Atom type</td><td>Partner ID</td><td>Partner x</td>"
		<< "<td> Partner y</td><td>Partner z</td><td>Distance</td></tr>\n";

	for(size_t i = 0; i < atoms.size(); i++)
	{
		for(size_t j = 0; j < atoms[i].bonds.size(); j++)
		{
			fout << "<tr><td>" << atoms[i].id << "</td><td>"
				<< atoms[i].type << "</td><td>"
				<< atoms[i].bonds[j].id << "</td><td>"
				<< atoms[i].bonds[j].coordinates.x << "</td><td>"
				<< atoms[i].bonds[j].coordinates.y << "</td><td>"
				<< atoms[i].bonds[j].coordinates.z << "</td><td>"
				<< atoms[i].bonds[j].dist << "</td></tr>\n";
		}
	}

	fout << "</table></body>\n";
	fout.close();
}

void Model::writeModel(string filename)
{
	ofstream fout(filename.c_str());

	fout << FILE_NAME << endl;
	fout << P_BONDS << " " << P_VASP << " " << WRITE_NN << endl;
	
	for(int k = 0; k < N_TYPES; k++){
		fout << nTypes[k] << " ";
	}
	fout << endl;
	
	for(int j = 0; j < N_TYPES; j++){
		fout << atomPercent[j] << " ";
	}
	fout << endl;
	
	for(int i = 0; i < N_TYPES; i++){
		fout << nAtoms[i] << " ";
	}
	fout << endl;
	
	fout << kT << endl;
	fout << I_NUM_SWITCHES << endl;
	fout << BOX_SIZE[0] << " " << BOX_SIZE[1] << " " << BOX_SIZE[2] << endl;
	fout << RELAX_V << " " << I_RELAX_V_SWITCHES << endl;

	for(size_t i = 0; i < atoms.size(); i++)
	{
		fout << atoms[i].id << "\t"
			<< atoms[i].type << "\t"
			<< atoms[i].coordinates.x << " "
			<< atoms[i].coordinates.y << " "
			<< atoms[i].coordinates.z << "\t"
			<< atoms[i].bonds.size() << "\t";
		for(size_t j = 0; j < atoms[i].bonds.size(); j++)
			fout << atoms[i].bonds[j].id << " "
			<< atoms[i].bonds[j].coordinates.x << " "
			<< atoms[i].bonds[j].coordinates.y << " "
			<< atoms[i].bonds[j].coordinates.z << "\t";
		fout << endl;
	}
	
	fout.close();
}

void Model::writeNearestNeighbors(string filename)
{
    ofstream fout(filename.c_str());

    fout << "<html><head><title>Atoms' Nearest Neighbors</title></head>\n"
         << "<body><table border=1>\n"
         << "<tr><td>Atom #</td><td>Nearest Neighbors</td>"
         << "<td>Distance</td></tr>\n";

    for(size_t j = 0; j < atoms.size(); j++)
    {
        fout << "<tr><td>" << j+1 << "</td>";

        // bonded neighbors
        for(size_t k = 0; k < atoms[j].bonds.size(); k++)
        {
            if(k == 0)
                fout << "<td> B " << atoms[j].bonds[k].id << "</td><td>"
                     << atoms[j].bonds[k].dist << "</td></tr>\n";
            else
                fout << "<tr><td></td><td> B " << atoms[j].bonds[k].id << "</td><td>"
                     << atoms[j].bonds[k].dist << "</td></tr>\n";
        }

        // non-bonded nearest neighbors
        for(size_t k = 0; k < atoms[j].closestNeighbors.size(); k++)
        {
            Point p = closestPointWithPBC(
                atoms[j].coordinates,
                atoms[atoms[j].closestNeighbors[k]-1].coordinates
            );

            double dist = p.distanceTo(atoms[j].coordinates);

            fout << "<tr><td></td><td>"
                 << atoms[atoms[j].closestNeighbors[k]-1].id
                 << "</td><td>" << dist << "</td></tr>\n";
        }
    }

    fout.close();
}


/* Sends the position data to a file
*/
void Model::writePositions(string filename)
{
	ofstream fout(filename.c_str());

	fout << "<html><head><title>Atom Positions</title></head>\n"
		<< "<body><table border=1>\n"
		<< "<tr><td>Atom ID</td><td>x</td><td>y</td><td>z</td></tr>\n";

	for(size_t i = 0; i < atoms.size(); i++)
	{
		fout << "<tr><td>" << atoms[i].id << "</td><td>"
			<< atoms[i].type << "</td><td>"
			<< atoms[i].coordinates.x << "</td><td>" 
			<< atoms[i].coordinates.y << "</td><td>" 
			<< atoms[i].coordinates.z << "</td></tr>\n";
	}

	fout << "</table></body>\n";
	fout.close();
}

/*Converts a model to .vasp format*/
void Model::writeVasp(string filename){

	keepInBox();

	ofstream fout(filename.c_str());

	fout << filename << endl << "1.00000000000000000" << endl;
	fout << BOX_SIZE[0] << " .0000000000000000 .0000000000000000" << endl;
	fout << ".0000000000000000 " << BOX_SIZE[1] << " .0000000000000000" << endl;
	fout << ".0000000000000000 .0000000000000000 " << BOX_SIZE[2] << endl;
	for(int i = 0; i < N_TYPES; i++){
		//	cout << nAtoms[i] << endl;
			fout << NAME[nTypes[i]] << " ";
	}
	fout << endl;
	for(int j = 0; j < N_TYPES; j++){

			fout << nAtoms[j] << " ";
	}
	fout << endl << "Direct" << endl;

	for(size_t i = 0; i < atoms.size(); i++)
	{
		fout << atoms[i].coordinates.x/BOX_SIZE[0] << " " 
			<< atoms[i].coordinates.y/BOX_SIZE[1] << " " 
			<< atoms[i].coordinates.z/BOX_SIZE[2] << "\n";
	}
	
	fout.close();

	return;
}

/* Sets the default values for constants
These can be edited at any time, but most should not
*/
void Model::setDefaultConstants()
{
	NAME = vector<string>(5, "Si");
	MASS = vector<double>(5, 1E-27);
	NUM_BONDS = vector<unsigned int>(5, 0);
	K_B = vector<vector<double> >(5, vector<double>(5, 0));
	B_0 = vector<vector<double> >(5, vector<double>(5, 0));
	D_0 = vector<vector<double> >(5, vector<double>(5, 0));
	COS_ANGLE = vector<double>(5, 0);			
	GAMMA = vector<vector<double> >(5, vector<double>(5, 0));
	K_OMEGA = vector<vector<vector<double> > >(5, vector<vector<double> >(5, vector<double>(5, 0)));

	SHELL_ONLY = false; 
	BOX_SIZE[0] = BOX_SIZE[1] = BOX_SIZE[2]  = 10.86; 
	SHELL_SIZE = 1.5;

	NUM_BONDS[0] = 4;
	NUM_BONDS[1] = 1;

	K_OMEGA[0][0][0] = 3.58;
	K_OMEGA[0][0][1] = 4.2;
	K_OMEGA[0][1][0] = 4.2;
	K_OMEGA[0][1][1] = 4.8;

	COS_ANGLE[0] = -1/3.0;

	K_B[0][0] = 9.08; 
	K_B[0][1] = 40;
	K_B[1][0] = 40;

	B_0[0][0] = 2.3513;
	B_0[0][1] = 1.5;
	B_0[1][0] = 1.5;

	D_0[0][0] = 2;
	D_0[0][1] = 2;
	D_0[1][0] = 2;
	D_0[1][1] = 2; 

	GAMMA[0][0] = 0.1;
	GAMMA[0][1] = 0; 
	GAMMA[1][0] = 0; 
	GAMMA[1][1] = 0.1;

	MASS[1] = 28 * 1.66E-27; 
	REPEL_DISTANCE = 4.2; 
}


