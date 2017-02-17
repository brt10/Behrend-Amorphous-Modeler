#Behrend Amorphous Modelling

###BAM:
The Behrend Amorphous Modeller takes an input file ending in `.input` and will generate output files with the following extensions:
* `.e` - energy file
* `.b` - bond file
* `.nn` - nearest neighbors file
* `.vasp` - VASP file
      
###Input File Format

```text
output_file_prefix=	Character	      Default= input_file_prefix	description= suffix for all ouput file names
energy=			T or F		Default= True			description= file with list of model energy for each step 
error=			T or F		Defualt= True			description= error messages 
output=			T or F		Default= True			description= file format with input file format for final model 
bond=	      		T or F		Default= True			description= file with atom and its number of bonds listed
vasp=       		T or F		Default= True			description= file with final model in vasp format
nn=		      	T or F		Default= True			description= file with atom and its neighbors listed	
atom_type_1=		Character	      Default= none			description= one or two letter character corresponding to atom from periodic table  	example= Si 
atom_type_2=		same as above
atom_number_1=		integer		Defualt= none			description= number of atoms of type one
atom_number_2= 		same as above
bond_switch_prob1=	real		Default= none		      	description= probability atom1 chosen for bond switch; value must be between 0 and 100
bond_switch_prob2=	same as above
atom_switch_code=       T or F		Default=False			description= function to switch atom types
atoms_to_switch=        two strings       Default= none			description= the first atom is switched for the second  example=  Si C 
atoms_switch_prob=	real		      Default=none			description= the probability that bond change is used instead of bond switch; must be between 0 and 100
temperature=		real	      	Default=none			description= temperature in KT units of eV
number_switches=  	integer		Defualt=0			      description= number of switching events created
lattice_constants=	real		      Default=0 0 0			description= three reals space seperated for x,y,z lattice lengths in Angstroms; 0 0 0 create lattice based on number and size of atoms
relax_volume=		T or F		Default=T		      	description= relax volume of model
volume_relax_time=	integer		Default=0	      		description= number of accepted switches before volume relaxed
atoms_fixed=		character	      Defualt=none			description= atom name to fix the position of atoms
lattice_constants=	real		      Default=0 0 0			description= three reals space seperated for x,y,z lattice lengths in Angstroms; 0 0 0 create lattice based on number and size of atoms
```
If lattice_constants given, then the format for atoms is:

`atom_type atom_number  x1 y1 z1  nbonds  IDa Xa Ya Za  IDb Xb Yb Zb  IDc Xc Yc Zc  IDd Xd Yd Zd`

From line 10 to end.

Line 9 contains, from left to right, the flag for volume relaxations, the number of accepted switches to wait before every volume relaxation, the first n atoms to fix the position (not moving), the first n atoms to totally exclude from swapping, the first n atoms to exclude from swapping each other (swap with other atoms is allowed).
Set the first flag to 1 to relax the volume of the current model and 0 to not relax the volume.  If this number is anything other than 0 or 1, an error message will be printed and the program will stop.  If the number of bond switches is not an integer or if this number is greater than the number of bond switches to be performed, an error message will be printed and the program will stop.



The output files will have the same basename as the input file.
