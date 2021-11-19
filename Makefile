MAIN=leviosam
CXX=g++
CXX_FLAGS=--std=c++11 -lpthread
LIB=-lsdsl -lhts -lz
CFLAGS=-g -Wall -O2
CLIB=-lz -lhts

OBJS = src/leviosam.o src/bam_aux.o src/bam_md.o src/bed.o src/chain.o src/collate.o src/leviosam_utils.o 
CDEPS= src/bam.h
DEPS = src/leviosam.hpp src/chain.hpp src/chain.hpp src/leviosam_utils.hpp src/collate.hpp src/robin_hood.h src/IITree.h
OLIB = liblvsam.a
TESTS = src/leviosam_test.cpp src/chain_test.cpp

all: $(MAIN) $(OLIB)

.c.o: %.c $(CDEPS)
	    $(CC) -c -o $@ $< $(CLIB) $(CFLAGS)

.cpp.o: %.cpp $(DEPS)
	    $(CXX) -c -o $@ $< $(LIB) $(CXX_FLAGS)

$(OLIB): $(OBJS)
		ar rcs $@ $^

$(MAIN): $(OBJS)
	    $(CXX) -o $@ $^ $(LIB) $(CXX_FLAGS)

# Tests
leviosam_test: src/leviosam_test.cpp
	    $(CXX) -o $@ $^ $(OLIB) $(LIB) $(CXX_FLAGS) -lgtest_main -lgtest

chain_test: src/chain_test.cpp
	    $(CXX) -o $@ $^ $(OLIB) $(LIB) $(CXX_FLAGS) -lgtest_main -lgtest

test: leviosam_test chain_test
		cd testdata; echo $(pwd); ../leviosam_test; ../chain_test; cd ../

# Clean
clean:
	    rm -f src/*.gch src/*.o $(MAIN) $(OLIB) leviosam_test chain_test
