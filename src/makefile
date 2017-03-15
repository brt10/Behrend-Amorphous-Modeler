CFLAGS = -Wall -O3 -std=c++11
LFLAGS = -Wall

bam: main.o model.o atom.o
	g++ $(LFLAGS) main.o -o bam
	ln -s ../src/bam ../bin/bam

main.o: main.cpp
	g++ $(CFLAGS) -o main.o -c main.cpp

model.o: model.cpp model.h
	g++ $(CFLAGS) -o model.o -c model.cpp model.h

atom.o: atom.cpp atom.h
	g++ $(CFLAGS) -o atom.o -c atom.cpp atom.h
