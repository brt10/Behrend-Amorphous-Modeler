/* Alicia Klinvex
April 22, 2008

Last modified by Brittany Hoard on June 16, 2012

Stores a set of atoms
Contains functions for amorphizing the atoms
*/

#include "Model.h"
//#include "Path.h"
#include <omp.h>
//#include <pthread.h>
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

inline double min(double a, double b)
{
	if(a < b) return a;
	else return b;
}

inline double max(double a, double b)
{
	if(a > b) return a;
	else return b;
}

/* Default constructor
*/
Model::Model()
{
	// seeds the random number generator
	srand((unsigned int)time(0));
	atoms.clear();
	IGNORE_ATOMS = 0; //completely ignore the first n atoms from swaps
	IGNORE_SWAP = 0;  //freeze the bonds between the first n atoms (bonds with other atoms can swap)
	FIX_ATOMS = 0;    //fix the first n atoms' position
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
	SHELL_SIZE = 4;
	framePause = 5;
	kT = 0.1;
	MAX_BONDS=4;
	I_NUM_SWITCHES = 0;
	I_RELAX_V_SWITCHES = 0;
	REPEL_DISTANCE = 7.2;
	FORCE_CONVERSION = 1.602E-29;
	CLUSTER_MAX_SHELL = 4;

	RADIAL_ON = true;
	ANGULAR_ON = true;
	REPEL_ON = true;

	//setDefaultConstants();
}

/* Sets the default values for constants
These can be edited at any time, but most should not
*/
// void Model::setDefaultConstants()
// {
// 	NAME = vector<string>(6, "Si");
// 	MASS = vector<double>(6, 1E-27);
// 	NUM_BONDS = vector<unsigned int>(6, 0);
// 	K_B = vector<vector<double> >(6, vector<double>(6, 0));
// 	B_0 = vector<vector<double> >(6, vector<double>(6, 0));
// 	D_0 = vector<vector<double> >(6, vector<double>(6, 0));
// 	COS_ANGLE = vector<double>(6, 0);
// 	GAMMA = vector<vector<double> >(6, vector<double>(6, 0));
// 	K_OMEGA = vector<vector<vector<double> > >(6, vector<vector<double> >(6, vector<double>(6, 0)));
//
// 	BOX_SIZE[0] = BOX_SIZE[1] = BOX_SIZE[2]  = 10.86;
//
// 	NUM_BONDS[0] = 4;
// 	NUM_BONDS[1] = 1;
//
// 	K_OMEGA[0][0][0] = 3.58;
// 	K_OMEGA[0][0][1] = 4.2;
// 	K_OMEGA[0][1][0] = 4.2;
// 	K_OMEGA[0][1][1] = 4.8;
//
// 	COS_ANGLE[0] = -1./3.0;
//
// 	K_B[0][0] = 9.08;
// 	K_B[0][1] = 40.;
// 	K_B[1][0] = 40.;
//
// 	B_0[0][0] = 3.04;
// 	B_0[0][1] = 3.04;
// 	B_0[1][0] = 3.04;
//
// 	D_0[0][0] = 2.;
// 	D_0[0][1] = 2.;
// 	D_0[1][0] = 2.;
// 	D_0[1][1] = 2.;
//
// 	GAMMA[0][0] = 0.1;
// 	GAMMA[0][1] = 0.;
// 	GAMMA[1][0] = 0.;
// 	GAMMA[1][1] = 0.1;
//
// 	MASS[1] = 28. * 1.66E-27;
//
// }

/* Performs numSwaps annealing transpositions
   Sends energy information to a file with a given filename, if one is specified
   Otherwise, there is no output, which is much faster
*/
void Model::performAnnealing(const int numSwaps, const int numVSwaps, string filename)
{
	bool OUTPUT_ON = true;
	ofstream fout;
	int acceptedSwaps = 0;

	if(OUTPUT_ON){
		fout.open(filename.c_str());
		fout << "#nswap initial_energy  relaxation_steps  final_energy  random_number  Metropolis_acceptance_probability\n";
	}

	double initEnergy, finEnergy, metProb, rand, minEnergy;
	vector<Atom> minAtoms = atoms, backup;
	int numSteps;

	char progress[4] = {'-','\\','|','/'};
	int progress_step=(numSwaps<400?1:int(numSwaps/400));

	if(REPEL_ON) initClosestNeighbors(); //initializes ClosestNeighbors for the first time
	minEnergy = initEnergy = getEnergy();//initial Energy

	if(numSwaps == 0){
		relax();
		if(RELAX_V == 1) relaxVolume();
		finEnergy = getEnergy();
		if(finEnergy < minEnergy){
			minEnergy = finEnergy;
			//	minAtoms = atoms;
		}
		if(OUTPUT_ON){
			fout << 0 << "  " << initEnergy << "  "
			<< numSteps << "  " << finEnergy << "  "
			<< rand << "  " << metProb << "\n";
		}
	}
	else for(int i = 0 ; i < numSwaps; i++)
	{
		if( i%progress_step == 0) cout << "\r" << progress[(i/progress_step)%4]
						<< " Progress: " << (i*100)/numSwaps << " %  "
						<< "   accepted swaps: " << acceptedSwaps ; fflush(stdout);

		backup = atoms;
		initEnergy = getEnergy();

		randomSwap();
		numSteps = relax();

		finEnergy = getEnergy();
		metProb = getMetropolisProb(initEnergy, finEnergy);
		rand = getRandomDecimal();

		if(rand >= metProb) atoms = backup;// rejected
		else{
			acceptedSwaps++;
			if(RELAX_V == 1 && acceptedSwaps%numVSwaps == 0){
				cout << "\n";
				relaxVolume();
				finEnergy = getEnergy();
			}
			if(finEnergy < minEnergy){ //If the (accepted) final energy is less than the current minimum, make this finEnergy the new minEnergy and preserve atom configuration
				minEnergy = finEnergy;
			//	minAtoms = atoms;
			}
			if(OUTPUT_ON){
				fout << i << "  " << initEnergy << "  "
					<< numSteps << "  " << finEnergy << "  "
					<< rand << "  " << metProb << "\n";
			}
		}

	}//end swaps

	cout << "\r  Progress: 100 %  " << "   accepted swaps: " << acceptedSwaps << endl ; fflush(stdout);

	fout << "#Total Number of accepted swaps: " << acceptedSwaps << " , minimum energy: " << minEnergy << "\n";
	fout.close();
}


/* Sets the forces on the atoms
*/
void Model::setForces()
{
	static size_t n = atoms.size(),i, j, k, jSubs, kSubs;
	static double multiplier, tempNum, cosijk, tempJ, tempK, dist;
	static Point pj,pk;

	// reset all forces
	for(i = 0; i < n; i++)
	{
// 		if((SHELL_ONLY && atoms[i].shellNumber > CLUSTER_MAX_SHELL))
// 			continue;

		atoms[i].force.x = 0;
		atoms[i].force.y = 0;
		atoms[i].force.z = 0;
	}

	#pragma omp parallel for default(shared) private(i,j,k,jSubs,kSubs,multiplier, tempNum, cosijk, tempJ, tempK, dist, pj, pk) schedule(static)
	for(i = 0; i < n; i++)
	{
// 		if((SHELL_ONLY && atoms[i].shellNumber > CLUSTER_MAX_SHELL))
// 			continue;

		for(j = 0; j < atoms[i].bonds.size(); j++)
		{
			jSubs = atoms[i].bonds[j].id;
			pj = closestPointWithPBC(atoms[i].coordinates, atoms[jSubs].coordinates);
			// set radial forces of two partner atoms at a time
			// twice as efficient, since it takes the "equal and opposite" rule into account
			if(RADIAL_ON  &&  i < jSubs)
			{
				multiplier = (K_B[atoms[i].type][atoms[jSubs].type] * (atoms[i].bonds[j].dist-B_0[atoms[i].type][atoms[jSubs].type]) / atoms[i].bonds[j].dist);

				tempNum = multiplier * (pj.x - atoms[i].coordinates.x);
				atoms[i].force.x += tempNum;
				#pragma omp atomic
				atoms[jSubs].force.x -= tempNum;

				tempNum = multiplier * (pj.y - atoms[i].coordinates.y);
				atoms[i].force.y += tempNum;
				#pragma omp atomic
				atoms[jSubs].force.y -= tempNum;

				tempNum = multiplier * (pj.z - atoms[i].coordinates.z); //was (atoms[i].bonds[j].coordinates.z - atoms[i].coordinates.z)
				atoms[i].force.z += tempNum;
				#pragma omp atomic
				atoms[jSubs].force.z -= tempNum;
			} // end radial force

			// set angular force of three atoms at a time
			// i is the "middle" atom and j/k are endpoints of the angle
			for(k = j+1; k < atoms[i].bonds.size() && ANGULAR_ON; k++)
			{
				kSubs = atoms[i].bonds[k].id;
				pk = closestPointWithPBC(atoms[i].coordinates, atoms[kSubs].coordinates);
				// avoids recomputing these repeatedly
				cosijk = ((pj.x - atoms[i].coordinates.x) * (pk.x - atoms[i].coordinates.x) +
					(pj.y - atoms[i].coordinates.y) * (pk.y - atoms[i].coordinates.y) +
					(pj.z - atoms[i].coordinates.z) * (pk.z - atoms[i].coordinates.z)) /
					(atoms[i].bonds[j].dist * atoms[i].bonds[k].dist);
				multiplier = K_OMEGA[atoms[i].type][atoms[jSubs].type][atoms[kSubs].type] *
					(COS_ANGLE[atoms[i].type]-cosijk);

				// find force on j
				tempJ = multiplier * ((pk.x-atoms[i].coordinates.x) / atoms[i].bonds[k].dist +
					(atoms[i].coordinates.x-pj.x) * cosijk / atoms[i].bonds[j].dist) / atoms[i].bonds[j].dist;
				#pragma omp atomic
				atoms[jSubs].force.x += tempJ;
				atoms[i].force.x -= tempJ;

				tempJ = multiplier * ((pk.y-atoms[i].coordinates.y) / atoms[i].bonds[k].dist +
					(atoms[i].coordinates.y-pj.y) * cosijk / atoms[i].bonds[j].dist) / atoms[i].bonds[j].dist;
				#pragma omp atomic
				atoms[jSubs].force.y += tempJ;
				atoms[i].force.y -= tempJ;

				tempJ = multiplier * ((pk.z-atoms[i].coordinates.z) / atoms[i].bonds[k].dist +
					(atoms[i].coordinates.z-pj.z) * cosijk / atoms[i].bonds[j].dist) / atoms[i].bonds[j].dist;
				#pragma omp atomic
				atoms[jSubs].force.z += tempJ;
				atoms[i].force.z -= tempJ;

				// find force on k
				tempK = multiplier * ((pj.x-atoms[i].coordinates.x) / atoms[i].bonds[j].dist +
					(atoms[i].coordinates.x-pk.x) * cosijk / atoms[i].bonds[k].dist) / atoms[i].bonds[k].dist;
				#pragma omp atomic
				atoms[kSubs].force.x += tempK;
				atoms[i].force.x -= tempK;

				tempK = multiplier * ((pj.y-atoms[i].coordinates.y) / atoms[i].bonds[j].dist +
					(atoms[i].coordinates.y-pk.y) * cosijk / atoms[i].bonds[k].dist) / atoms[i].bonds[k].dist;
				#pragma omp atomic
				atoms[kSubs].force.y += tempK;
				atoms[i].force.y -= tempK;

				tempK = multiplier * ((pj.z-atoms[i].coordinates.z) / atoms[i].bonds[j].dist +
					(atoms[i].coordinates.z-pk.z) * cosijk / atoms[i].bonds[k].dist) / atoms[i].bonds[k].dist;
				#pragma omp atomic
				atoms[kSubs].force.z += tempK;
				atoms[i].force.z -= tempK;
			} // end angular force
		} // end j loop

		if(REPEL_ON)
		{
			// repulsive force term
			for(j = 0; j < atoms[i].closestNeighbors.size(); j++)
			{
//				jSubs = atoms[atoms[i].closestNeighbors[j]].id;
				// find closest version of atom
				pj = closestPointWithPBC(atoms[i].coordinates, atoms[atoms[i].closestNeighbors[j]].coordinates);
				dist = pj.distanceTo(atoms[i].coordinates);

				// if the distance is greater than D_0, there is no force
				if(dist < D_0[atoms[i].type][atoms[atoms[i].closestNeighbors[j]].type])
				{
					multiplier = 3 * GAMMA[atoms[i].type][atoms[atoms[i].closestNeighbors[j]].type] * pow(D_0[atoms[i].type][atoms[atoms[i].closestNeighbors[j]].type] - dist, 2) / dist;
					tempJ=multiplier * (atoms[i].coordinates.x - pj.x);
					atoms[i].force.x += tempJ;
					#pragma omp atomic
					atoms[jSubs].force.x -= tempJ;  //TEST MYMOD: no force was applied to jsubs

					tempJ = multiplier * (atoms[i].coordinates.y - pj.y);
					atoms[i].force.y += tempJ;
					#pragma omp atomic
					atoms[jSubs].force.y -= tempJ;

					tempJ = multiplier * (atoms[i].coordinates.z - pj.z);
					atoms[i].force.z += tempJ;
					#pragma omp atomic
					atoms[jSubs].force.z -= tempJ;
				}
			}
		} // end if repel on
	} // end main loop
}

/* Sets the radial energy of the atoms
*/
void Model::setRadialEnergy()
{
	/* radial energy is divided between the two atoms */
	static size_t i,j,n=atoms.size(),jSubs;

	#pragma omp parallel for default(shared) private(i,j,jSubs) schedule(static)
	for(i = 0; i < n; i++)
	{

// 		if(SHELL_ONLY && atoms[i].shellNumber > CLUSTER_MAX_SHELL+1)
// 			continue;

		atoms[i].radialEnergy = 0;
		for(j = 0; j < atoms[i].bonds.size(); j++)
		{
			jSubs = atoms[i].bonds[j].id;
			// divided by 4 not 2 because two atoms share the energy
			atoms[i].radialEnergy += K_B[atoms[i].type][atoms[jSubs].type] * pow(atoms[i].bonds[j].dist - B_0[atoms[i].type][atoms[jSubs].type], 2) / 4;
		}
	}
}

/* Sets the angular energy of the system (angular energy belongs to the center atom)
*/
void Model::setAngularEnergy()
{
	static size_t i,j,k,n=atoms.size(),jSubs,kSubs;
	static Point pj,pk;

	// select atom
	#pragma omp parallel for default(shared) private(i,j,k,jSubs,kSubs,pj,pk) schedule(static)
	for(i = 0; i < n; i++)
	{
// 		if(SHELL_ONLY && atoms[i].shellNumber > CLUSTER_MAX_SHELL+1)
// 			continue;

		atoms[i].angularEnergy = 0;
		for(j = 0; j < atoms[i].bonds.size(); j++) // select first bond
		{
			jSubs = atoms[i].bonds[j].id;
			pj = closestPointWithPBC(atoms[i].coordinates, atoms[jSubs].coordinates);
			for(k = j+1; k < atoms[i].bonds.size(); k++) // select second bond
			{
				kSubs = atoms[i].bonds[k].id;
				pk = closestPointWithPBC(atoms[i].coordinates, atoms[kSubs].coordinates);
				// 0.5 0.647 2.35^2 (cosijk + 1/3)^2
				atoms[i].angularEnergy += (K_OMEGA[atoms[i].type][atoms[jSubs].type][atoms[kSubs].type] *
					pow(cosijk(atoms[i].coordinates, pj, pk, atoms[i].bonds[j].dist, atoms[i].bonds[k].dist) - COS_ANGLE[atoms[i].type], 2) / 2);
			}
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

	#pragma omp parallel for default(shared) private(i,j,tempPt,tempDist) schedule(static)
	for(i = 0; i < n; i++)
	{
// 		if(SHELL_ONLY && atoms[i].shellNumber > CLUSTER_MAX_SHELL+1)
// 			continue;

		atoms[i].repulsiveEnergy = 0;
		for(j = 0; j < atoms[i].closestNeighbors.size(); j++)
		{
			// find closest version of this point that might be repelling
			tempPt = closestPointWithPBC(atoms[i].coordinates, atoms[atoms[i].closestNeighbors[j]].coordinates);
			tempDist = tempPt.distanceTo(atoms[i].coordinates);

			if(tempDist < D_0[atoms[i].type][atoms[atoms[i].closestNeighbors[j]].type]) // add term to energy ~ 0.5 (3.84 - |rm - rn|)^3
				atoms[i].repulsiveEnergy += (GAMMA[atoms[i].type][atoms[atoms[i].closestNeighbors[j]].type] * pow(D_0[atoms[i].type][atoms[atoms[i].closestNeighbors[j]].type] - tempDist, 3) / 2);
		}
	}
}


/* Allows the atoms to move until the forces become very small
*/
int Model::relax()
{
	static int numSteps;
	static double biggestForce;

	for(numSteps = 0; ; numSteps++)
	{
		setForces();
		setAccelerations();
		biggestForce = getLargestForce();

		if( biggestForce < 0.1 ) break;
		else if ( numSteps >= 300)
			{/*cout << "\nbiggestForce = " << biggestForce << " , numSteps = " << numSteps << endl;*/break;}
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

		keepInBox();
	}

	return numSteps;
}


/*Change the box size of an existing model, preserving the atoms' fractional coordinates
*/
void Model::changeBoxSize(double variation){

	static size_t i,j,n=atoms.size();

	#pragma omp parallel for default(shared) private(i,j) schedule(static)
	for(i = 0; i < n; i++)
	{
		atoms[i].coordinates.x *= variation;
		atoms[i].coordinates.y *= variation;
		atoms[i].coordinates.z *= variation;

		for(j = 0; j < atoms[i].bonds.size(); j++){
			atoms[i].bonds[j].coordinates.x *= variation;
			atoms[i].bonds[j].coordinates.y *= variation;
			atoms[i].bonds[j].coordinates.z *= variation;
		}

	}

	//Change box size
	BOX_SIZE[0] *= variation;
	BOX_SIZE[1] *= variation;
	BOX_SIZE[2] *= variation;

//	if(REPEL_ON) initClosestNeighbors(); //MYMOD: commented out because ClosestNeighbors don't depend on the cell size
	for(int i = 0; i < 50 && relax() > 0; i++); //relax atoms

	return;
}


void Model::relaxVolume(){ //MYMOD: now can deal with orthorombic cells

	static double scaleup = 1.01;
	static double scaledown = 1./scaleup;
	int nchange=0;

	double origSize0 = BOX_SIZE[0];
	double origSize1 = BOX_SIZE[1];
	double origSize2 = BOX_SIZE[2];

	double origEnergy = getEnergy();
	double newEnergy;

	//try to increase the volume
	cout << "\rchangeBoxSize: " << ++nchange << "   " ; fflush(stdout);
	changeBoxSize(scaleup);
	newEnergy = getEnergy();

	if(newEnergy < origEnergy){ //keep increasing the volume

		origEnergy = newEnergy;
		origSize0 = BOX_SIZE[0];
		origSize1 = BOX_SIZE[1];
		origSize2 = BOX_SIZE[2];

		while(1){

			cout << "\rchangeBoxSize: " << ++nchange << "   " ; fflush(stdout);
			changeBoxSize(scaleup);
			newEnergy = getEnergy();

			if(newEnergy < origEnergy){
				origEnergy = newEnergy;
				origSize0 = BOX_SIZE[0];
				origSize1 = BOX_SIZE[1];
				origSize2 = BOX_SIZE[2];
			}
			else{ //step back and exit
				cout << "\rchangeBoxSize: " << --nchange << "   " ; fflush(stdout);
				BOX_SIZE[0] = origSize0;
				BOX_SIZE[1] = origSize1;
				BOX_SIZE[2] = origSize2;
				break;
			}
		}
	}

	//try to reduce the volume
	else if(newEnergy > origEnergy){

		//restore the box size
		cout << "\rchangeBoxSize: " << --nchange << "   " ; fflush(stdout);
		BOX_SIZE[0] = origSize0;
		BOX_SIZE[1] = origSize1;
		BOX_SIZE[2] = origSize2;

		//and reduce the volume
		cout << "\rchangeBoxSize: " << --nchange << "   " ; fflush(stdout);
		changeBoxSize(scaledown);
		newEnergy = getEnergy();

		if(newEnergy < origEnergy){

			//keep the new box
			origEnergy = newEnergy;
			origSize0 = BOX_SIZE[0];
			origSize1 = BOX_SIZE[1];
			origSize2 = BOX_SIZE[2];

			//keep reducing until energy decrease
			while(1){

				cout << "\rchangeBoxSize: " << --nchange << "   " ; fflush(stdout);
				changeBoxSize(scaledown);
				newEnergy = getEnergy();

				if(newEnergy < origEnergy){
					origEnergy = newEnergy;
					origSize0 = BOX_SIZE[0];
					origSize1 = BOX_SIZE[1];
					origSize2 = BOX_SIZE[2];
				}
				else{
					cout << "\rchangeBoxSize: " << ++nchange << "   " ; fflush(stdout);
					BOX_SIZE[0] = origSize0;
					BOX_SIZE[1] = origSize1;
					BOX_SIZE[2] = origSize2;
					break;
				}
			}
		}
		else { //restoring the initial box size
			cout << "\rchangeBoxSize: " << ++nchange << "   " ; fflush(stdout);
			BOX_SIZE[0] = origSize0;
			BOX_SIZE[1] = origSize1;
			BOX_SIZE[2] = origSize2;
		}
	}

	else{ //restoring the initial box size
		cout << "\rchangeBoxSize: " << --nchange << "   " ; fflush(stdout);
		BOX_SIZE[0] = origSize0;
		BOX_SIZE[1] = origSize1;
		BOX_SIZE[2] = origSize2;
	}

	cout << "\n";

	return;

}

/* Sets the accelerations based on the forces
*/
void Model::setAccelerations()
{
	static size_t n = atoms.size(),i;

	for(i = FIX_ATOMS; i < n; i++)
	{
		// converts the force to kg A/fs^2, then divides by the atom mass in kg
		atoms[i].acceleration.x = atoms[i].force.x * FORCE_CONVERSION / MASS[atoms[i].type];
		atoms[i].acceleration.y = atoms[i].force.y * FORCE_CONVERSION / MASS[atoms[i].type];
		atoms[i].acceleration.z = atoms[i].force.z * FORCE_CONVERSION / MASS[atoms[i].type];
	}
}


/* Sets the positions of the atoms based on the acceleration
*/
void Model::setPositions(const double timestep)
{
	static Vector displacement;
	double multiplier = pow(timestep, 2) / 2.0;
	static size_t n = atoms.size(),i,j,k,jSubs;
	static Point pj;

	// updates the position of non-fixed atoms
	#pragma omp parallel for default(shared) private(i,displacement) schedule(static)
	for(i = FIX_ATOMS; i < n; i++)
	{
// 		if((SHELL_ONLY && atoms[i].shellNumber > CLUSTER_MAX_SHELL) )
// 			continue;

		// calculate displacement
		displacement.x = atoms[i].acceleration.x * multiplier;
		displacement.y = atoms[i].acceleration.y * multiplier;
		displacement.z = atoms[i].acceleration.z * multiplier;
		// adjust coordinates
		atoms[i].coordinates.x += displacement.x;
		atoms[i].coordinates.y += displacement.y;
		atoms[i].coordinates.z += displacement.z;
	}

	// updates ALL the bonds
	#pragma omp parallel for default(shared) private(i,j,k,jSubs,displacement,pj) schedule(static)
	for(i = 0; i < n; i++)
	{
		for(j = 0; j < atoms[i].bonds.size(); j++)
		{
			// adjust bond length i-j
			jSubs = atoms[i].bonds[j].id;
			pj = closestPointWithPBC(atoms[i].coordinates, atoms[jSubs].coordinates);
			atoms[i].bonds[j].dist = atoms[i].coordinates.distanceTo(pj);

			// adjust bond length j-i
			for(k = 0; k < atoms[jSubs].bonds.size(); k++)
			{
				if(atoms[i].id == atoms[jSubs].bonds[k].id) //OMP atomization should not be required here because bond j-i changes only when i changes
				{
					atoms[jSubs].bonds[k].coordinates.x = atoms[i].coordinates.x; //MYMOD: to take into account FIX_ATOMS
					atoms[jSubs].bonds[k].coordinates.y = atoms[i].coordinates.y;
					atoms[jSubs].bonds[k].coordinates.z = atoms[i].coordinates.z;
					atoms[jSubs].bonds[k].dist = atoms[i].bonds[j].dist;
				}
			}
		}
	}
}


/* Performs a random WWW bond transposition on the system
*/
void Model::randomSwap()//WARNING: when IGNORE_ATOMS is not small this routine becomes dramatically inefficient!
{
	static size_t  aSubs, bSubs, cSubs, dSubs, initR, trialR, n=atoms.size(), randIndex;
	static unsigned int randNum;
	static unsigned int maxPercent;
	static size_t i,j,*numtype,**typelist;
	static bool once=true;

	//enumerates the atoms for every type. This is executed only at the first call.
	if(once){
		typelist = new size_t*[N_TYPES]; for(i = 0; i < N_TYPES; i++) typelist[i] = new size_t[n]; // typelist[NTYPES][n]
		numtype = new size_t[N_TYPES]; // numtype[NTYPES]

		//fills the array
		for(i = 0; i < N_TYPES; i++){
			numtype[i]=0;
			for(j = IGNORE_ATOMS; j < n; j++){
				if(atoms[j].type == nTypes[i]){
					typelist[i][numtype[i]]=atoms[j].id;
					numtype[i]++;
				}
			}
		}
		once=false;
	}

	//cout << "\nSelecting bonds to swap... " ; fflush(stdout);

	while(true)
	{
		// choose A based on atom swap probabilities
		maxPercent = 0;
		if(N_TYPES > 1)
		{
			randNum = getRandomInteger(); //gets a random integer [1..100]
			for(i = 0; i < N_TYPES; i++)
			{
				maxPercent = maxPercent + atomPercent[i];
				if(randNum <= maxPercent)
				{ // type selected; selecting one atom of type i and exit loop
					randIndex = ( rand() % numtype[i] );
					aSubs = typelist[i][randIndex];
					break;
				}
			}
		}
		else
		{
			aSubs = (rand() % numtype[0]);
		}

		// choose B
		initR = rand() % atoms[aSubs].bonds.size();
		bSubs = atoms[aSubs].bonds[initR].id;
		if( bSubs < IGNORE_ATOMS ) continue;
		
		// choose C
		initR = rand() % atoms[bSubs].bonds.size();
		cSubs = atoms[bSubs].bonds[initR].id;

		if(cSubs == aSubs || isNearestNeighbor(cSubs, atoms[aSubs]))
		{
			// so try another C
			trialR = (initR) % atoms[bSubs].bonds.size();
			cSubs = atoms[bSubs].bonds[trialR].id;

			while((cSubs == aSubs || isNearestNeighbor(cSubs, atoms[aSubs])) && trialR != initR)
			{
				trialR = (trialR) % atoms[bSubs].bonds.size();
				cSubs = atoms[bSubs].bonds[trialR].id;
			}

			// no value of C will work
			if(trialR == initR)
				continue;
		}
		if( cSubs < IGNORE_ATOMS ) continue;
		// b and c can be both < IGNORE_SWAP, since the b-c bond is preserved
		// a can be < IGNORE_SWAP only if both b,c are not
		if( aSubs < IGNORE_SWAP && (bSubs < IGNORE_SWAP || cSubs < IGNORE_SWAP) )
			continue;

		// choose D
		initR = rand() % atoms[cSubs].bonds.size();
		dSubs = atoms[cSubs].bonds[initR].id;

		// there is a problem with D
		if(dSubs == bSubs || isNearestNeighbor(dSubs, atoms[bSubs]))
		{
			// so try another D
			trialR = (initR) % atoms[cSubs].bonds.size();
			dSubs = atoms[cSubs].bonds[trialR].id;

			while((dSubs == bSubs || isNearestNeighbor(dSubs, atoms[bSubs])) && trialR != initR)
			{
				trialR = (trialR) % atoms[cSubs].bonds.size();
				dSubs = atoms[cSubs].bonds[trialR].id;
			}

			// no value of D will work
			if(trialR == initR)
				continue;
		}
		if( dSubs < IGNORE_ATOMS ) continue;

		// b and c can be both < IGNORE_SWAP, since the b-c bond is preserved
		// d can be < IGNORE_SWAP only if both b,c are not
		if( dSubs < IGNORE_SWAP && (bSubs < IGNORE_SWAP || cSubs < IGNORE_SWAP) )
			continue;

		break;
	}

	// swapping them to be a-c-b-d
	//cout << "Swapping " << aSubs+1 << " " << bSubs+1 << " " << cSubs+1 << " " << dSubs+1 << endl; fflush(stdout);
	swapBonds(aSubs, bSubs, cSubs, dSubs);
}


/* Swaps the bonds of four atoms based on the WWW algorithm
*/
inline void Model::swapBonds(const size_t aSubs, const size_t bSubs, const size_t cSubs, const size_t dSubs)
{
	static size_t i;
	for(i = 0; i < MAX_BONDS; i++)
	{
		Point temp;
		// a-b becomes a-c
		if(i < atoms[aSubs].bonds.size() && atoms[aSubs].bonds[i].id == bSubs)
		{
			temp = closestPointWithPBC(atoms[aSubs].coordinates, atoms[cSubs].coordinates);
			atoms[aSubs].bonds[i] = Partner(temp, cSubs, temp.distanceTo(atoms[aSubs].coordinates));
		}

		// b-a becomes b-d
		if(i < atoms[bSubs].bonds.size() && atoms[bSubs].bonds[i].id == aSubs)
		{
			temp = closestPointWithPBC(atoms[bSubs].coordinates, atoms[dSubs].coordinates);
			atoms[bSubs].bonds[i] = Partner(temp, dSubs, temp.distanceTo(atoms[bSubs].coordinates));
		}

		// c-d becomes c-a
		if(i < atoms[cSubs].bonds.size() && atoms[cSubs].bonds[i].id == dSubs)
		{
			temp = closestPointWithPBC(atoms[cSubs].coordinates, atoms[aSubs].coordinates);
			atoms[cSubs].bonds[i] = Partner(temp, aSubs, temp.distanceTo(atoms[cSubs].coordinates));
		}

		// d-c becomes d-b
		if(i < atoms[dSubs].bonds.size() && atoms[dSubs].bonds[i].id == cSubs)
		{
			temp = closestPointWithPBC(atoms[dSubs].coordinates, atoms[bSubs].coordinates);
			atoms[dSubs].bonds[i] = Partner(temp, bSubs, temp.distanceTo(atoms[dSubs].coordinates));
		}

	}

	if(REPEL_ON) initClosestNeighbors(); //updates ClosestNeighbors
	//if(SHELL_ONLY) initShellNumbers(aSubs, bSubs, cSubs, dSubs);
}


/* Sends the bond information to a file
*/
void Model::writeBonds(string filename)
{
	size_t i,j,n=atoms.size();
	ofstream fout(filename.c_str());

	fout << "#Atom_ID Atom_type Partner_ID Partner_x Partner_y Partner_z Distance\n";

	for(i = 0; i < n; i++)
	{
		for(j = 0; j < atoms[i].bonds.size(); j++)
		{
			fout << atoms[i].id << " "
				<< atoms[i].type << " "
				<< atoms[i].bonds[j].id << " "
				<< atoms[i].bonds[j].coordinates.x << " "
				<< atoms[i].bonds[j].coordinates.y << " "
				<< atoms[i].bonds[j].coordinates.z << " "
				<< atoms[i].bonds[j].dist << "\n";
		}
	}

	fout.close();
}

/*  writes the .anneal file
 */
void Model::writeModel(string filename)
{
	ofstream fout(filename.c_str());
	size_t n=atoms.size();

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
	fout << RELAX_V << " " << I_RELAX_V_SWITCHES << " " << FIX_ATOMS << " " << IGNORE_ATOMS << " " << IGNORE_SWAP << endl;

	for(size_t i = 0; i < n; i++)
	{
		fout << atoms[i].id+1 << "\t"
			<< atoms[i].type << "\t"
			<< atoms[i].coordinates.x << " "
			<< atoms[i].coordinates.y << " "
			<< atoms[i].coordinates.z << "\t"
			<< atoms[i].bonds.size() << "\t";
		for(size_t j = 0; j < atoms[i].bonds.size(); j++)
			fout << atoms[i].bonds[j].id+1 << " "
			<< atoms[i].bonds[j].coordinates.x << " "
			<< atoms[i].bonds[j].coordinates.y << " "
			<< atoms[i].bonds[j].coordinates.z << "\t";
		fout << endl;
	}

	fout.close();
}

/*  writes the .nn file
 */
void Model::writeNearestNeighbors(string filename)
{
	ofstream fout(filename.c_str());
	size_t n=atoms.size(),j,k;
	double dist;
	Point p;

	fout << "<html><head><title>Atoms' Nearest Neighbors</title></head>\n"
		 << "<body><table border=1>\n"
		 << "<tr><td>Atom #</td><td>Nearest Neighbors</td>"
		 << "<td>Distance</td></tr>\n";

	for(j = 0; j < n; j++){
		fout << "<tr><td>" << j+1 << "</td>";
		for(k = 0; k < atoms[j].bonds.size(); k++){
			if(k == 0)
				fout << "<td>" << " B " << atoms[j].bonds[k].id+1 << "</td><td>" << atoms[j].bonds[k].dist << "</td></tr>\n";
			else
				fout << "<tr><td><td>" << " B " << atoms[j].bonds[k].id+1 << "</td><td>" << atoms[j].bonds[k].dist << "</td></tr>\n";
		}

		for(k = 0; REPEL_ON && k < atoms[j].closestNeighbors.size(); k++){
			p = closestPointWithPBC(atoms[j].coordinates, atoms[atoms[j].closestNeighbors[k]].coordinates);
			dist = p.distanceTo(atoms[j].coordinates);
			fout << "<tr><td><td>" << atoms[atoms[j].closestNeighbors[k]].id+1 << "</td><td>" << dist << "</td></tr>\n";
		}

	}

	fout.close();

}

/* Sends the position data to a file
*/
void Model::writePositions(string filename)
{
	size_t n=atoms.size();
	ofstream fout(filename.c_str());

	fout << "<html><head><title>Atom Positions</title></head>\n"
		<< "<body><table border=1>\n"
		<< "<tr><td>Atom ID</td><td>x</td><td>y</td><td>z</td></tr>\n";

	for(size_t i = 0; i < n; i++)
	{
		fout << "<tr><td>" << atoms[i].id+1 << "</td><td>"
			<< atoms[i].type << "</td><td>"
			<< atoms[i].coordinates.x << "</td><td>"
			<< atoms[i].coordinates.y << "</td><td>"
			<< atoms[i].coordinates.z << "</td></tr>\n";
	}

	fout << "</table></body>\n";
	fout.close();
}

/*   writes the vasp file
 */
void Model::writeVasp(string filename){

	size_t n=atoms.size();
	keepInBox();

	ofstream fout(filename.c_str());

	fout << filename << endl << "1.0000000000000000" << endl;
	fout << BOX_SIZE[0] << " .0000000000000000 .0000000000000000" << endl;
	fout << ".0000000000000000 " << BOX_SIZE[1] << " .0000000000000000" << endl;
	fout << ".0000000000000000 .0000000000000000 " << BOX_SIZE[2] << endl;
	for(int i = 0; i < N_TYPES; i++){
			fout << NAME[nTypes[i]] << " ";
	}
	fout << endl;
	for(int j = 0; j < N_TYPES; j++){

			fout << nAtoms[j] << " ";
	}
	fout << endl << "Direct" << endl;

	for(size_t i = 0; i < n; i++)
	{
		fout << atoms[i].coordinates.x/BOX_SIZE[0] << " "
			<< atoms[i].coordinates.y/BOX_SIZE[1] << " "
			<< atoms[i].coordinates.z/BOX_SIZE[2] << "\n";
	}

	fout.close();

	return;
}

/* Finds the version of atom p2 that's closest to atom p1
*/
inline Point Model::closestPointWithPBC(const Point &p1, const Point &p2)
{
	static Point bestPoint;
// 	int dx, dy, dz;
// 	int	startX = int(p2.x / BOX_SIZE[0]) - 1,
// 		startY = int(p2.y / BOX_SIZE[1]) - 1,
// 		startZ = int(p2.z / BOX_SIZE[2]) - 1;
// 	double minDistance;
// 	Point temp;
//
// 	if(p2.x < 0) startX--;
// 	if(p2.y < 0) startY--;
// 	if(p2.z < 0) startZ--;
//
// 	for(dx = startX, minDistance=INT_MAX; dx <= 1; dx++)
// 	{
// 		temp.x = p2.x + dx * BOX_SIZE[0]; //x box size
//
// 		if(fabs(temp.x - p1.x) < minDistance)
// 		{
// 			minDistance = fabs(temp.x - p1.x);
// 			bestPoint.x = temp.x;
// 		}
// 	}
//
// 	for(dy = startY, minDistance=INT_MAX; dy <= 1; dy++)
// 	{
// 		temp.y = p2.y + dy * BOX_SIZE[1];//y box size
//
// 		if(fabs(temp.y - p1.y) < minDistance)
// 		{
// 			minDistance = fabs(temp.y - p1.y);
// 			bestPoint.y = temp.y;
// 		}
// 	}
//
// 	for(dz = startZ, minDistance=INT_MAX; dz <= 1; dz++)
// 	{
// 		temp.z = p2.z + dz * BOX_SIZE[2];//z box size
//
// 		if(fabs(temp.z - p1.z) < minDistance)
// 		{
// 			minDistance = fabs(temp.z - p1.z);
// 			bestPoint.z = temp.z;
// 		}
// 	}

	bestPoint.x = p2.x - round( (p2.x-p1.x) / BOX_SIZE[0] )*BOX_SIZE[0];
	bestPoint.y = p2.y - round( (p2.y-p1.y) / BOX_SIZE[1] )*BOX_SIZE[1];
	bestPoint.z = p2.z - round( (p2.z-p1.z) / BOX_SIZE[2] )*BOX_SIZE[2];

	return bestPoint;
}

/* computes the cosine of an angle using the dot product
*/
inline double Model::cosijk(const Point &i, const Point &j, const Point &k, const double rij, const double rik)
{
	return (((j.x - i.x) * (k.x - i.x) +
		(j.y - i.y) * (k.y - i.y) +
		(j.z - i.z) * (k.z - i.z)) /
		(rij * rik));
}

/* Sets and returns the energy of the system
*/
inline double Model::getEnergy()
{
	static size_t n = atoms.size(),i;
	static double rEnergy, angEnergy, rpEnergy;

	rEnergy = angEnergy = rpEnergy = 0.;

	if(RADIAL_ON)	setRadialEnergy();
	if(ANGULAR_ON)	setAngularEnergy();
	if(REPEL_ON)	setRepulsiveEnergy();

	for(i = 0; i < n; i++){
		if(RADIAL_ON) rEnergy = rEnergy + atoms[i].radialEnergy;
		if(ANGULAR_ON) angEnergy = angEnergy + atoms[i].angularEnergy;
		if(REPEL_ON) rpEnergy = rpEnergy + atoms[i].repulsiveEnergy;
	}

	return (rEnergy+angEnergy+rpEnergy);
}

/* Returns the largest force on an atom in the system
*/
inline double Model::getLargestForce()
{
	static size_t n = atoms.size(),i;
	static double maxForce;
	maxForce = 0.;

	for(i = FIX_ATOMS; i < n; i++)
	{
// 		if((SHELL_ONLY && atoms[i].shellNumber > CLUSTER_MAX_SHELL))
// 			continue;

		// using absolute value
		maxForce = max(maxForce, fabs(atoms[i].force.x));
		maxForce = max(maxForce, fabs(atoms[i].force.y));
		maxForce = max(maxForce, fabs(atoms[i].force.z));
	}

	return maxForce;
}

/* Returns the metropolis acceptance probability
*/
inline double Model::getMetropolisProb(const double Einitial, const double Efinal)
{
	if(kT == 0.) return double(Einitial>=Efinal);  // if kT==0 the argument of the exponential gives infinity or -infinity,
	else return min(1, exp((Einitial-Efinal)/kT)); // the exponential gives infinity or zero, and the function returns 1 or 0
}

/* Returns a random decimal in the range (0, 1)
*/
inline double Model::getRandomDecimal()
{
	return (((rand()%9999)+1)/10000.0);
}

/* Returns a random integer in the range (1, 100)
*/
inline unsigned int Model::getRandomInteger()
{
	return ((rand()%100)+1);
}

/* Initializes the atoms nearby contributing to the repulsive energy
*/
void Model::initClosestNeighbors()
{
	static size_t n=atoms.size(),i,j;
	double tempDist;

	for(i = 0; i < n; i++)
	{
// 		if(SHELL_ONLY && atoms[i].shellNumber > CLUSTER_MAX_SHELL+1)
// 			continue;

		// reset array of nearby atoms
		atoms[i].closestNeighbors.clear();

		for(j = 0; j < n; j++)
		{
// 			if(SHELL_ONLY && atoms[j].shellNumber > CLUSTER_MAX_SHELL+1)
// 				continue;

			// ignore the atom if it is a neighbor of the current one
			if(i == j || isNearestNeighbor(j, atoms[i]) || isSecondNearestNeighbor(j, i)) continue;

			tempDist = atoms[i].coordinates.distanceTo(closestPointWithPBC(atoms[i].coordinates, atoms[j].coordinates));

			if(tempDist <= REPEL_DISTANCE) atoms[i].closestNeighbors.push_back(j);
		}
	}
}

/* Returns true if atoms[id] is bonded to atom
*/
inline bool Model::isNearestNeighbor(const size_t id, const Atom &atom)
{
	static size_t i;

	for(i = 0; i < atom.bonds.size(); i++)
	{
		if(atom.bonds[i].id == id)
			return true;
	}
	return false;
}

/* Returns true if there exists an atom p such that
   atoms[id] is bonded to p and atoms[currentSubs] is bonded to p
*/
inline bool Model::isSecondNearestNeighbor(const size_t id, const size_t currentSubs)
{
	static size_t neighborID,i;

	for(i = 0; i < atoms[currentSubs].bonds.size(); i++)
	{
		neighborID = atoms[currentSubs].bonds[i].id;

		// if id is a nearest neighbor to one of current's partners,
		// id and current are secondNearestNeighbors
		if(isNearestNeighbor(id, atoms[neighborID]))
			return true;
	}
	return false;
}

/* Puts the atoms back in the periodic box
*/
inline void Model::keepInBox()
{
	static size_t i,j,n=atoms.size();

	for(i = 0; i < n; i++)
	{
		while(atoms[i].coordinates.x < 0)
		{
				atoms[i].coordinates.x += BOX_SIZE[0];
				for(j = 0; j < atoms[i].bonds.size(); j++)
					atoms[i].bonds[j].coordinates.x += BOX_SIZE[0];
		}
		while(atoms[i].coordinates.x >= BOX_SIZE[0])
		{
				atoms[i].coordinates.x -= BOX_SIZE[0];
				for(j = 0; j < atoms[i].bonds.size(); j++)
					atoms[i].bonds[j].coordinates.x -= BOX_SIZE[0];
		}

		while(atoms[i].coordinates.y < 0)
		{
			atoms[i].coordinates.y += BOX_SIZE[1];
			for(j = 0; j < atoms[i].bonds.size(); j++)
				atoms[i].bonds[j].coordinates.y += BOX_SIZE[1];
		}
		while(atoms[i].coordinates.y >= BOX_SIZE[1])
		{
			atoms[i].coordinates.y -= BOX_SIZE[1];
			for(j = 0; j < atoms[i].bonds.size(); j++)
				atoms[i].bonds[j].coordinates.y -= BOX_SIZE[1];
		}

		while(atoms[i].coordinates.z < 0)
		{
			atoms[i].coordinates.z += BOX_SIZE[2];
			for(j = 0; j < atoms[i].bonds.size(); j++)
				atoms[i].bonds[j].coordinates.z += BOX_SIZE[2];
		}
		while(atoms[i].coordinates.z >= BOX_SIZE[2])
		{
			atoms[i].coordinates.z -= BOX_SIZE[2];
			for(j = 0; j < atoms[i].bonds.size(); j++)
				atoms[i].bonds[j].coordinates.z -= BOX_SIZE[2];
		}
	}
}

/* Reads Constants From a Model File
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
	fin >> RELAX_V >> RELAX_V_SWITCHES >> FIX_ATOMS >> IGNORE_ATOMS >> IGNORE_SWAP;


	if(P_BONDS != 0 && P_BONDS != 1){
		cerr << "The leftmost number on line 2 (.b file creation) in the file " << C_FILE_NAME << " should be 0 or 1.  Program aborted." << endl;
		error = true;
	}

	if(P_VASP != 0 && P_VASP != 1){
		cerr << "The middle number on line 2 (.vasp file creation) in the file " << C_FILE_NAME << " should be 0 or 1.  Program aborted." << endl;
		error = true;
	}

	if(WRITE_NN != 0 && WRITE_NN != 1){
		cerr << "The rightmost number on line 2 (.nn file creation) in the file " << C_FILE_NAME << " should be 0 or 1.  Program aborted." << endl;
		error = true;
	}

	for(int j = 0; j < N_TYPES; j++){
		if(nTypes[j] >= NUM_ELEMENTS){
			cerr << "You have entered an invalid atom type on line 3.  Program aborted." << endl;
			error = true;
		}

	}

	if(N_TYPES == 1 && nTypes[0] == 2){
		cerr << "You have entered only H atoms on line 3.  You cannot have only H atoms in a model.  Program aborted." << endl;
		error = true;
	}

	for(int j = 0; j < N_TYPES; j++){
		if(nAtoms[j] == 0){
			cerr << "You entered 0 atoms for one or more atom types on line 5.  Please remove the corresponding atom type(s) from line 3.  Program aborted." << endl;
			error = true;
		}
	}

	unsigned int totalPercent = 0;
	for(int j = 0; j < atomPercent.size(); j++){
		totalPercent = totalPercent + atomPercent[j];
	}

	if(totalPercent != 100){
		cerr << "The numbers on line 4 should add up to 100.  Program aborted." << endl;
		error = true;
	}

	if(N_TYPES != nAtoms.size()){
		cerr << "The number of atom types on line 3 does not match the numbers of atoms on line 5.  Program aborted." << endl;
		error = true;
	}

	if(kT < 0 || kT > 20){
		cerr << "The kT value (line 6) in the file " << C_FILE_NAME << " needs to be between 0 and 20.  Program aborted." << endl;
		error = true;
	}

	if (static_cast<int>(NUM_SWITCHES) == NUM_SWITCHES){
		I_NUM_SWITCHES = static_cast<int>(NUM_SWITCHES);
	}
	else{
		cerr << "The number of bond swaps (line 7) in the file " << C_FILE_NAME << " needs to be an integer value.  Program aborted." << endl;
		error = true;
	}

	if(RELAX_V != 0 && RELAX_V != 1){
		cerr << "The leftmost number on line 9 in the file " << C_FILE_NAME << " should be 0 or 1.  Program aborted." << endl;
		error = true;
	}

	if (static_cast<int>(RELAX_V_SWITCHES) == RELAX_V_SWITCHES){
		I_RELAX_V_SWITCHES = static_cast<int>(RELAX_V_SWITCHES);
	}
	else{
		cerr << "The number of bond swaps before volume relaxation (line 9) in the file " << C_FILE_NAME << " needs to be an integer value.  Program aborted." << endl;
		error = true;
	}

	if(RELAX_V != 0 && I_RELAX_V_SWITCHES > I_NUM_SWITCHES){
		cerr << "The number of bond swaps before volume relaxation (line 9) in the file " << C_FILE_NAME << " cannot be greater than the total number of bond swaps performed.  Program aborted." << endl;
		error = true;
	}

	int id, numBonds, type, id2;
	double x, y, z, x2, y2, z2;
	MAX_BONDS=0;
	//Reads a model from a file. WARNING: Expecting atom id from 1 to N in input, while stores atom id from 0 to N-1 in memory.
	if(!error){
		while(fin >> id >> type >> x >> y >> z >> numBonds){

			if(numBonds > MAX_BONDS) MAX_BONDS = numBonds;
			atoms.push_back(Atom(Point(x, y, z), id-1));
			atoms[atoms.size()-1].type = type;

			for(int i = 0; i < numBonds; i++)
			{
				fin >> id2 >> x2 >> y2 >> z2;
				atoms[atoms.size()-1].bonds.push_back(Partner(Point(x2,y2,z2), id2-1, atoms[atoms.size()-1].coordinates.distanceTo(closestPointWithPBC(Point(x,y,z),Point(x2,y2,z2)))));
			}
		}
	}
	fin.close();

	if(FIX_ATOMS>0 || IGNORE_ATOMS>0) cout << "WARNING: fixing the position (bonds) of the first " << FIX_ATOMS << " (" << IGNORE_ATOMS << ") atoms." << endl;

	if(error) return false;
	else return true;
}
