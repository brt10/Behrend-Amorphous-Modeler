/* Created by Alicia Klinvex
Modified by Joseph Synowka

Stores a set of atoms
Contains functions for amorphizing the atoms
*/
#ifndef MODEL_H
#define MODEL_H

//#include "Python.h"
#include "Point.h"
#include "Atom.h"
#include <string>

struct Model
{
	//constructor
	Model();
	bool bondAtoms();

	void changeBoxSize(double newBoxSize);

	//getters
	inline double getEnergy();
	inline double getLargestForce();

	//Create models
	void makeCrystal(const int numBoxes, const int e, const bool interface);
	void makeSICrystal(const int numBoxes, const unsigned int e);
	void insertH(const Point point, const int siId);
	bool makePreheatedCrystal(const int numAtoms);
	bool makePreheatedCrystal(const int numAtoms, const double x, const double y, const double z);
	bool makeDefectedSystem();
	bool randomPlacement(const int numAtoms, const int type, const double x, const double y, const double z, const double x_length, const double y_length, const double z_length);

	void deleteAtom(const int id);
	void deleteElement(const int type);
	bool makeInterface(string baseFile, const int type, const int numAtoms, bool fixed, bool ignore);
	bool bondWithInterface(int bondTo);

	//simulation functions
	inline double getMetropolisProb(const double Einitial, const double Efinal);
	inline double getRandomDecimal();
	inline unsigned int getRandomInteger();
	inline double cosijk(const Point &i, const Point &j, const Point &k, const double rij, const double rik);
	inline Point closestPointWithPBC(const Point &p1, const Point &p2);

	void performAnnealing(const int numSwaps, const int numVSwaps, string filename);
	void randomSwap();
	void randomSwap2();
	int relax();
	void relaxVolume();
	inline void swapBonds(const size_t aSubs, const size_t bSubs, const size_t cSubs, const size_t dSubs);
	void initAtomsInShell();
	void initClosestNeighbors();

	inline bool isNearestNeighbor(const size_t id, const Atom &atom);
	inline bool isSecondNearestNeighbor(const size_t id, const size_t currentSubs);
	inline void keepInBox();
	void initShellNumbers(const size_t aSubs, const size_t bSubs, const size_t cSubs, const size_t dSubs);
	void fourNearestNeighbors(const char *filename);
	void keepInInterfaceBox(double oldX);

	//writing functions
	void writeBonds(string filename);
	void writeModel(string filename);
	void writeNearestNeighbors(string filename);
	void writeNewModel(string filename);
	void writePositions(string filename);
	void writeVasp(string filename);

	//reading functions
	void readConstants(const char* filename);
	bool readInput();
	void readModel();

	//setters
	void setAngularEnergy();
	void setAccelerations();
	void setEnergy();
	void setForces();
	void setPositions(const double timestep);
	void setRadialEnergy();
	void setRepulsiveEnergy();
	void setDefaultConstants();

	//data members
	vector<Atom> atoms;
	vector<unsigned int> nTypes;
	vector<unsigned int> nAtoms;
	vector<unsigned int> nIds;
	vector<unsigned int> atomPercent;
	unsigned int NEW_ATOMS;

	unsigned int NUM_ELEMENTS;
	unsigned int MAX_BONDS;
	double kT;
	float NUM_SWITCHES;
	unsigned int I_NUM_SWITCHES;
	float RELAX_V_SWITCHES;
	unsigned int I_RELAX_V_SWITCHES;
	bool SHELL_ONLY;
	bool IGNORE_SI_SWAP;
	unsigned int IGNORE_ATOMS;
	unsigned int IGNORE_SWAP;
	unsigned int FIX_ATOMS;
	unsigned int P_BONDS;
	unsigned int P_VASP;
	unsigned int WRITE_NN;
	unsigned int N_TYPES;
	unsigned int CREATE_MODEL;
	unsigned int SI_3;
	unsigned int SI_4;
	unsigned int SI_5;
	unsigned int NUM_H;
	unsigned int RELAX_V;
	string FILE_NAME;
	string C_FILE_NAME;
	string OUTFILE_NAME;

	double BOX_SIZE[3]; //[0]= x length, [1]= y length, [2]= z length
	double SHELL_SIZE;

	vector<string> NAME;
	vector<unsigned int> NUM_BONDS;
	vector<vector<vector<double> > > K_OMEGA;
	vector<double> COS_ANGLE;
	vector<vector<double> > K_B;
	vector<vector<double> > B_0;
	vector<vector<double> > D_0;
	vector<vector<double> > GAMMA;
	vector<double> MASS;

	double REPEL_DISTANCE;
	int CLUSTER_MAX_SHELL;
	double FORCE_CONVERSION;
	double framePause;
	double OLD_X;

	bool RADIAL_ON;
	bool ANGULAR_ON;
	bool REPEL_ON;
	bool INTER;
};

#endif
