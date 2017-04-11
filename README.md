# Behrend Amorphous Modelling

### BAM:
The Behrend Amorphous Modeller takes an input file ending in `.input` and will generate output files with the same prefix and the following extensions:
* `.e` - energy file
* `.b` - bond file
* `.nn` - nearest neighbors file
* `.vasp` - VASP file
      
### Input File Format

```
output_file_prefix=	string				Default= input file prefix	description= prefix for all ouput file names
energy=			true or false			Default= true			description= file with list of model energy for each step 
error=			true or false			Defualt= true			description= error messages 
output=			true or false			Default= true			description= file format with input file format for final model 
bond=	      		true or false			Default= true			description= file with atom and its number of bonds listed
vasp=       		true or false			Default= true			description= file with final model in vasp format
nn=		      	true or false			Default= true			description= file with atom and its neighbors listed	
atom=			string, int, double		Default= none			description= string for the atom type, int for the quantity, and double for the switching probability  	example= Si 45 0.37
atom= 			same as above
atom=			same as above
...
atom_switch_code=       true or false			Default= false			description= function to switch atom types
atoms_to_switch=        two strings     		Default= none			description= the first atom is switched for the second  example=  Si C 
atoms_switch_prob=	double				Default= none			description= the probability that bond change is used instead of bond switch; must be between 0 and 100
temperature=		double	  		    	Default= none			description= temperature in KT units of eV
number_switches=  	integer				Defualt= 0			description= number of switching events created
lattice_constants=	double				Default= 0 0 0			description= three reals space seperated for x,y,z lattice lengths in Angstroms; 0 0 0 create lattice based on number and size of atoms
relax_volume=		true or false			Default= true		      	description= relax volume of model
volume_relax_time=	integer				Default= 0	      		description= number of accepted switches before volume relaxed
atoms_fixed=		string				Defualt= none			description= atom name to fix the position of atoms
```
If lattice_constants given, then the format for atoms is:

`atom_type atom_number  x1 y1 z1  nbonds  IDa Xa Ya Za  IDb Xb Yb Zb  IDc Xc Yc Zc  IDd Xd Yd Zd`

