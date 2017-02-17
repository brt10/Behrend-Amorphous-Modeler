CFLAGS = -Wall -O3
LFLAGS = -Wall

bam: main.o
	g++ $(LFLAGS) main.o -o bam

main.o: main.cpp
	g++ $(CFLAGS) -o main.o -c main.cpp