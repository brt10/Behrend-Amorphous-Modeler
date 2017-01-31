#include <omp.h>
#include "Model.h"
#include <time.h>
#include <iostream>
#include <fstream>
#include <string>
using namespace std;

// This is main
int main(int argc, char * argv[])
{
#ifdef OMP
	{	int nthreads, tid;
		#pragma omp parallel private(tid)
		{
			tid=omp_get_thread_num();
			if (tid == 0) {
				nthreads = omp_get_num_threads();
				printf("Number of threads = %d\n", nthreads);
			}
		}
        }
#else
	printf("Serial version\n");
#endif

   Model m; //creates a default model with default constants
   if(argc > 1){
		string str(argv[1]);
		m.C_FILE_NAME = str;
   }
   else { 
	   cout << "No input file specified.  Specify file." << endl;
	   cin >> m.C_FILE_NAME;
   }

   m.readConstants("Tuttle_constants.txt");
   if(m.readInput()){

//		string startFileName = m.FILE_NAME  + ".start";

// 		if(m.atoms.size() == 0){
// 			m.C_FILE_NAME = startFileName;
// 			m.makeDefectedSystem();
// 			m.writeModel(startFileName);
//
// 		}

		cout << "checkpoint 1";

		if(m.I_NUM_SWITCHES >= 0){
			string energyFileName = m.FILE_NAME + ".e";
			string annealedFileName = m.FILE_NAME + ".anneal";
			m.C_FILE_NAME = energyFileName;
			m.performAnnealing(m.I_NUM_SWITCHES, m.I_RELAX_V_SWITCHES, energyFileName);
			m.C_FILE_NAME = annealedFileName;
			m.writeModel(annealedFileName);
/*			atoms = minAtoms;
                        string annealedFileName = m.FILE_NAME + ".minenergy";
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

// 		if(m.RELAX_V == 1 && m.I_NUM_SWITCHES == 0){
// 			m.relaxVolume();
// 			m.writeModel(startFileName);
// 		}

   }

   return 0;
}
