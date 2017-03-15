#include "model.h"
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <map>
#include <vector>

using namespace std;

struct atom {
	string atom_type;
	int atom_number, nbonds, IDa, IDb, IDc, IDd;
	double x1, y1, z1, Xa, Ya, Za, Xb, Yb, Zb, Xc, Yc, Zc, Xd, Yd, Zd;
};

int main(int argc, char* argv[]) {

	if (2 != argc) {
		cerr << "USAGE: bam file.input" << endl;
		return 1;
	}

	string configLine; // string buffer for reading in the config file
	bool eOutput, errOutput, mOutput, bOutput, vOutput, nnOutput; // Options for choosing which output files to generate
	bool switchCode, relaxVolume;
	string atomTypeOne, atomTypeTwo, atomFixed, atomToSwitchOne, atomToSwitchTwo; // strings for atoms, could be 2 char arrays
	int atomQuantityOne, atomQuantityTwo, numberSwitches, relaxTime;
	double bondSwitchProbOne, bondSwitchProbTwo, atomSwitchProb, temp, latticeConsts[3];
	string basename; // name of output files

	ifstream configFile(argv[1]);

	string key, value, buffer;
	map<string, string> configMap;
	vector<atom> configAtoms;
	configAtoms.reserve(100);
	while (getline(configFile, buffer)) {
		if ('#' == buffer[0])
			continue;
		else if (!buffer.find('=')) {
			atom* n = new atom;
			configFile >> n->atom_type >> n->atom_number >> n->x1 >> n->y1 >> n->z1 >> n->nbonds >> n->IDa >> n->Xa >> n->Ya >> \
n->Za >> n->IDb >> n->Xb >> n->Yb >> n->Zb >> n->IDc >> n->Xc >> n->Yc >> n->Zc >> n->IDd >> n->Xd >> n->Yd >> n->Zd;
			configAtoms.push_back(*n);
		}
		else {
			key = buffer.substr(0, buffer.find('='));
			value = buffer.erase(0, buffer.find('=') + 1);
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
	istringstream(configMap["atom_type_1"]) >> atomTypeOne;
	istringstream(configMap["atom_type_2"]) >> atomTypeTwo;
	istringstream(configMap["atom_number_1"]) >> atomQuantityOne;
	istringstream(configMap["atom_number_2"]) >> atomQuantityTwo;
	istringstream(configMap["bond_switch_prob1"]) >> bondSwitchProbOne;
	istringstream(configMap["bond_switch_prob2"]) >> bondSwitchProbTwo;
	istringstream(configMap["atom_switch_code"]) >> boolalpha >> switchCode;
	istringstream(configMap["atoms_to_switch"]) >> atomToSwitchOne >> atomToSwitchTwo;
	istringstream(configMap["atoms_switch_prob"]) >> atomSwitchProb;
	istringstream(configMap["temperature"]) >> temp;
	istringstream(configMap["number_switches"]) >> numberSwitches;
	istringstream(configMap["lattice_constants"]) >> latticeConsts[0] >> latticeConsts[1] >> latticeConsts[2];
	istringstream(configMap["relax_volume"]) >> boolalpha >> relaxVolume;
	istringstream(configMap["volume_relax_time"]) >> relaxTime;
	istringstream(configMap["atoms_fixed"]) >> atomFixed;

	cout << "basename is " << basename << endl;
	cout << "eOutput is " << eOutput << endl;
	cout << "errOutput is " << errOutput << endl;
	cout << "mOutput is " << mOutput << endl;
	cout << "bOutput is " << bOutput << endl;
	cout << "vOutput is " << vOutput << endl;
	cout << "nnOutput is " << nnOutput << endl;
	cout << "atomTypeOne is " << atomTypeOne << endl;
	cout << "atomTypeTwo is " << atomTypeTwo << endl;
	cout << "atomQuantityOne is " << atomQuantityOne << endl;
	cout << "atomQuantityTwo is " << atomQuantityTwo << endl;
	cout << "bondSwitchProbOne is " << bondSwitchProbOne << endl;
	cout << "bondSwitchProbTwo is " << bondSwitchProbTwo << endl;
	cout << "switchCode is " << switchCode << endl;
	cout << "atomToSwitchOne and Two are " << atomToSwitchOne << " " << atomToSwitchTwo << endl;
	cout << "atomSwitchProb is " << atomSwitchProb << endl;
	cout << "temp is " << temp << endl;
	cout << "numberSwitches is " << numberSwitches << endl;
	cout << "latticeConsts are " << latticeConsts[0] << " " << latticeConsts[1] << " " << latticeConsts[2] << endl;
	cout << "relaxVolume is " << relaxVolume << endl;
	cout << "relaxTime is " << relaxTime << endl;
	cout << "atomFixed is " << atomFixed << endl;

	return 0;
}
