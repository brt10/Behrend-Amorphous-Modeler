default:
	g++ -O3 -Wno-pragmas Model.cpp main.cpp -o bam.x
parallel:
	g++ -O3 -DOMP -fopenmp Model.cpp main.cpp -o bam.x 
serial:
	g++ -O3 -Wno-pragmas Model.cpp main.cpp -o bam.x

