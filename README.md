#Behrend Amorphous Modelling
#### bam

####You can edit the README with Markdown 
Check this for the basics in markdown https://help.github.com/articles/github-flavored-markdown/

###Todo:
    Move Todo list to seperate TODO file?
    To improve compatability/usability:
      -change file extensions to more standard formats (.txt or .csv possibly)
      -seperate files an end user would modify and input and output files from the actual src files
      -configure makefile to place bam.exe outside of src folder
      
###Input File Format

BAM takes input in the form of a file. The simplest input would be along the lines of:
```text
outputfilesprefix
1 1 1           //[bond data] [vasp file] [nn table -need to remove]
3               //# of atoms [type #1] [type #2] [etc. -check?] -located in tuttle_constants
100             //[% of bond switching -atom type 1] [% type 2]
50              //[number of atoms of type 1] [type 2] etc.
0.1             //[Boltsmann Cnst*Temperature (eV)]
1               //# bond switches           
0 0 0           //lattice constants: if 0 then auto]
1 1000 0 0      //[relax volume] [initial bond switches before vol. relax] [fixed atom (0 or N): fixes first N rows of atoms]                   // [no bond switching. (0,N)]
1 0 0 0 0          //input lines: 0 gives random atoms in box
                //otherwise # type x y z #nn #n1  --change so the code references the index numbers of the neighbors
```

This would be most likely inside a file such as **testinput.input**.
Note: This is described in https://github.com/awadYehya/bam/blob/master/src/My_input_manual.txt
Do We Need this?

###Output

