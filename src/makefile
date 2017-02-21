CFLAGS = -Wall -O3 -std=c++11
LFLAGS = -Wall

bam: main.o
	g++ $(LFLAGS) main.o -o bam

main.o: main.cpp
	g++ $(CFLAGS) -o main.o -c main.cpp
