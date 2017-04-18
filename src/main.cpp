#include "model.h"
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <map>
#include <tuple>
#include <vector>

using namespace std;

int main(int argc, char* argv[]) {

	if (2 != argc) {
		cerr << "USAGE: bam file.input" << endl;
		return 1;
	}

	string configLine; // string buffer for reading in the config file
	bool eOutput, errOutput, mOutput, bOutput, vOutput, nnOutput; // Options for choosing which output files to generate
	bool switchCode, relaxVolume;
	string atomTypeOne, atomTypeTwo, atomFixed, atomToSwitchOne, atomToSwitchTwo; // strings for atoms, could be 2 char arrays
	int numberSwitches, relaxTime;
	double atomSwitchProb, temp, latticeConsts[3];
	string basename; // name of output files

	ifstream configFile(argv[1]);

	string key, value, buffer;
	map<string, string> configMap; // dictionary for configuration file options
	map<string, tuple<int, double> > uAtoms; // u(nspecified)Atoms
	vector<atom> configAtoms; // config(uration file)Atoms
	configAtoms.reserve(100);
	while (getline(configFile, buffer) && !buffer.empty()) {
		if ('#' == buffer[0]) // if the first character of the buffer is a #, ignore it
			continue;
		else if (string::npos == buffer.find('=')) { // if '=' is not found in string buffer
			cout << "found an atom line: " << buffer << endl;
			istringstream stream(buffer);
			while (stream) {
				int tempID;       // intermediary variables for extracting data from string
				string tempType;
				double tempX, tempY, tempZ; //
				stream >> tempX >> tempY >> tempZ >> tempID >> tempType;
				atom *a = new atom(tempX, tempY, tempZ, tempID, tempType);
				configAtoms.push_back(*a);
				string leftover;
				stream >> leftover; // any unexpected string after the config variables
				if(!leftover.empty())
				    cerr << "\"" << leftover << "\" at the end of the atom line. Revise config file" << endl;
			}
		}
		else {
			key = buffer.substr(0, buffer.find('=')); // "key" will be from beginning of string to '='
			value = buffer.erase(0, buffer.find('=') + 1); // "value" will be everything after, so erase '=' and before
			if("atom" == key) {
			    string name;
			    int quantity;
			    double probability;
			    istringstream(value) >> name >> quantity >> probability; // read name, quantity and prob from atom= line
			    uAtoms[name] = tuple<int, double>(quantity, probability); // store the atom= line in the uAtoms map
			} else
			    configMap[key] = value;
		}
	}

	istringstream(configMap["output_file_prefix"]) >> basename;
	istringstream(configMap["energy"]) >> boolalpha >> eOutput;
	istringstream(configMap["error"]) >> boolalpha >> errOutput;
	istringstream(configMap["output"]) >> boolalpha >> mOutput;
	istringstream(configMap["bond"]) >> boolalpha >> bOutput;
	istringstream(configMap["vasp"]) >> boolalpha >> vOutput;
	istringstream(configMap["nn"]) >> boolalpha >> nnOutput;
	istringstream(configMap["atom_switch_code"]) >> boolalpha >> switchCode;
	istringstream(configMap["atoms_to_switch"]) >> atomToSwitchOne >> atomToSwitchTwo;
	istringstream(configMap["atoms_switch_prob"]) >> atomSwitchProb;
	istringstream(configMap["temperature"]) >> temp;
	istringstream(configMap["number_switches"]) >> numberSwitches;
	istringstream(configMap["lattice_constants"]) >> latticeConsts[0] >> latticeConsts[1] >> latticeConsts[2];
	istringstream(configMap["relax_volume"]) >> boolalpha >> relaxVolume;
	istringstream(configMap["volume_relax_time"]) >> relaxTime;
	istringstream(configMap["atoms_fixed"]) >> atomFixed;

	cout << "Configuration file options: " << endl;
	cout << "basename is " << basename << endl;
	cout << "eOutput is " << eOutput << endl;
	cout << "errOutput is " << errOutput << endl;
	cout << "mOutput is " << mOutput << endl;
	cout << "bOutput is " << bOutput << endl;
	cout << "vOutput is " << vOutput << endl;
	cout << "nnOutput is " << nnOutput << endl;
       	cout << "switchCode is " << switchCode << endl;
	cout << "atomToSwitchOne and Two are " << atomToSwitchOne << " " << atomToSwitchTwo << endl;
	cout << "atomSwitchProb is " << atomSwitchProb << endl;
	cout << "temp is " << temp << endl;
	cout << "numberSwitches is " << numberSwitches << endl;
	cout << "latticeConsts are " << latticeConsts[0] << " " << latticeConsts[1] << " " << latticeConsts[2] << endl;
	cout << "relaxVolume is " << relaxVolume << endl;
	cout << "relaxTime is " << relaxTime << endl;
	cout << "atomFixed is " << atomFixed << endl;

	cout << endl << "uAtoms: " << endl;
	for(const auto &m : uAtoms) {
	    cout << m.first << " " << get<0>(m.second) << " " << get<1>(m.second) << endl;
	}

	cout << endl << "configAtoms: " << endl;
	for (atom& a : configAtoms) {
		cout << a.getID() << " " << a.getType() << ": " << a.getX() << ", " << a.getY() << ", " << a.getZ() << endl;
	}

	model model(uAtoms, configAtoms, basename);
	
	return 0;
}
