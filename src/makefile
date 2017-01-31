CFLAGS = -Wall -O3
LFLAGS = -Wall
DEFAULTC = -DOMP -fopenmp
SERIALC = -Wno-pragmas

bam: main.o Model.o
	g++ $(LFLAGS) Model.o main.o -o bam

serial: main.serial.o Model.serial.o
	g++ $(LFLAGS) Model.serial.o main.serial.o -o serial

Model.o: Model.cpp Model.h
	g++ $(CFLAGS) $(DEFAULTC) -o Model.o -c Model.cpp

main.o: Model.cpp Model.h
	g++ $(CFLAGS) $(DEFAULTC) -o main.o -c main.cpp

Model.serial.o: Model.cpp Model.h
	g++ $(CFLAGS) $(SERIALC) -o Model.serial.o -c Model.cpp

main.serial.o: Model.cpp Model.h
	g++ $(CFLAGS) $(SERIALC) -o main.serial.o -c main.cpp

Model.cpp: Point.h Atom.h

Atom.cpp: Partner.h Vector.h
