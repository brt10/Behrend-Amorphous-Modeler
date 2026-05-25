/*#include <QApplication>
#include "modelwindow.h"

int main(int argc, char *argv[])
{
    QApplication app(argc, argv);
    ModelWindow window;
    window.show();
    return app.exec();
}*/
#include "Model.h"
#include <time.h>
#include <iostream>
#include <fstream>
#include <string>
using namespace std;

int main(int argc, char * argv[])
{

/**********TEST TO INSURE PROPER FUNCTION OF CODE**********/ 
   cout << "test1";
   Model m;
   if(argc > 1){
		string str(argv[1]);
		m.C_FILE_NAME = str;
   }


   m.readConstants("Tuttle_constants.txt");

   if(m.readInput()){
		string startFileName = m.FILE_NAME  + ".start";
		if(m.atoms.size() == 0){
			m.C_FILE_NAME = startFileName;
			m.makeDefectedSystem();
			m.writeModel(startFileName);
		}
		
		// SLT - Calculate initial energies and forces	
		string initialEnergyFileName = m.FILE_NAME + ".inite";
		string initialForcesFileName = m.FILE_NAME + ".initf";	
		m.writeInitialEnergy(initialEnergyFileName);
		m.writeInitialForces(initialForcesFileName);
		m.getEnergy();
		
		if(m.I_NUM_SWITCHES > 0){
			string energyFileName = m.FILE_NAME + ".e";
			string annealedFileName = m.FILE_NAME + ".anneal";
			m.C_FILE_NAME = energyFileName;
			m.performAnnealing(m.I_NUM_SWITCHES, m.I_RELAX_V_SWITCHES, energyFileName);
			m.C_FILE_NAME = annealedFileName;
			m.writeModel(annealedFileName);
/*			atoms = minAtoms;
                        string annealedFileName = m.FILE_NAME + ".minenergy";
                        m.C_FILE_NAME = annealedFileName;
                        m.writeModel(annealedFileName);
*/		}
		
		if(m.P_BONDS == 1){
			string bondsFileName = m.FILE_NAME + ".b";
			m.writeBonds(bondsFileName);
		}
		
		if(m.P_VASP == 1){
			string vaspFileName = m.FILE_NAME + ".vasp";
			m.writeVasp(vaspFileName);
		}

		if(m.WRITE_NN == 1){
			string nnFileName = m.FILE_NAME + ".nn";
			m.writeNearestNeighbors(nnFileName);
		}

		if(m.RELAX_V == 1 && m.I_NUM_SWITCHES == 0){
			m.relaxVolume();
			m.writeModel(startFileName);
		}
	}
   
/********************************************************/


   return 0;  
}
