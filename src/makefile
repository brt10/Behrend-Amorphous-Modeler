CFLAGS = -Wall -O3 -std=c++11
LFLAGS = -Wall

all: bam clean

bam: main.o model.o atom.o
	g++ $(LFLAGS) main.o model.o atom.o -o bam
	if [ ! -h ../bin/bam ]; then\
		ln -sf ../src/bam ../bin/bam;\
	fi

main.o: main.cpp model.h atom.h
	g++ $(CFLAGS) -o main.o -c main.cpp

model.o: model.cpp model.h
	g++ $(CFLAGS) -o model.o -c model.cpp

atom.o: atom.cpp atom.h
	g++ $(CFLAGS) -o atom.o -c atom.cpp

clean :
	rm -f main.o model.o atom.o
