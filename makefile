CPP=g++
INCLUDE=./include
CPPFLAG=-std=c++11 -I$(INCLUDE)
BINARY=bin

all: bin_dir pileup2bed

bin_dir:
	mkdir -p $(BINARY)

pileup2bed:
	$(CPP) -o bin/pileup2bed $(CPPFLAG) src/pileup2bed.cpp
