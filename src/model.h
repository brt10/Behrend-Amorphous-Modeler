#ifndef MODEL_H
#define MODEL_H

#include "atom.h"
#include <string>
#include <vector>

using namespace std;

class model {
	private:
		double boxSize[3]; // dimensions of the model box
		vector<atom> atoms;
	public:
		model(); // (default) constructor
		double getRandomDecimal();
		bool makeDefectedSystem();
		bool bondAtoms();
};

#endif
