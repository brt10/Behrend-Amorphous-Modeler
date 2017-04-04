#ifndef MODEL_H
#define MODEL_H

#include "atom.h"
#include <string>
#include <vector>

using namespace std;

class model {
	private:
		double boxSize[3]; // dimensions of the model box
		vector<atom> atoms; // vector of atoms in the model
	public:
		model(vector<int> uAtoms, vector<atom> configAtoms, string efn); // Constructor for random atoms + explicit atoms
		model(vector<int> uAtoms, string efn); // Constructor for random atoms
		double getRandomDecimal(); // returns a random decimal for randomly generating points
		bool bondAtoms(); // checks if the atoms bond properly
		void setBoxSize(double size); // sets the starting size of the box according to the number of atoms
		double getBoxSize(); // returns box size
};

#endif
