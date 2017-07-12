SHELL = /bin/sh
CC = g++
CFLAGS =  -static-libgcc -L. -w -Wall -Wextra -pedantic -Wshadow -funroll-loops -pthread -std=c++11 -DNDEBUG -O3 -march=native

all: aligner indexer main
aligner: sureMap.cpp
	$(CC) $(CFLAGS) -o SureMap-Aligner.o sureMap.cpp  -fopenmp
indexer: index_builder.cpp
	$(CC) $(CFLAGS) -o SureMap-Indexer.o ./psascan_src/utils.cpp index_builder.cpp  -ldivsufsort -ldivsufsort64 -fopenmp
clean:
	/bin/rm -f *.o
nuclear:
	/bin/rm -f psascan *.o
