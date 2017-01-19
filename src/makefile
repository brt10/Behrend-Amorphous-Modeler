CFLAGS = -Wall -O3
LFLAGS = -Wall
DEFAULTC = -DOMP -fopenmp
SERIALC = -Wno-pragmas

bam: main.o Model.o
	g++ $(LFLAGS) Model.o main.o -o bam

serial: main.serial.o Model.serial.o
	g++ $(LFLAGS) Model.serial.o main.serial.o -o bam

Model.o: Model.h
	g++ $(CFLAGS) $(DEFAULTC) -o Model.o -c Model.cpp

main.o: Model.h
	g++ $(CFLAGS) $(DEFAULTC) -o main.o -c main.cpp

Model.serial.o: Model.h
	g++ $(CFLAGS) $(SERIALC) -o Model.serial.o -c Model.cpp

main.serial.o: Model.h
	g++ $(CFLAGS) $(SERIALC) -o main.serial.o -c main.cpp

Model.h: Point.h Atom.h

Atom.h: Partner.h Vector.h
