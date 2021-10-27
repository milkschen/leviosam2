MAIN=leviosam
CXX=g++
CXX_FLAGS=--std=c++11 -lpthread
LIB=-lsdsl -lhts -lz
CFLAGS=-g -Wall -O2
CLIB=-lz -lhts

OBJS = src/leviosam.o src/chain.o src/leviosam_utils.o src/bam_aux.o src/bam_md.o
CDEPS= src/bam.h
DEPS = src/leviosam.hpp src/chain.hpp

all: $(MAIN)

.c.o: %.c $(CDEPS)
	    $(CC) -c -o $@ $< $(CLIB) $(CFLAGS)

.cpp.o: %.cpp $(DEPS)
	    $(CXX) -c -o $@ $< $(LIB) $(CXX_FLAGS)

$(MAIN): $(OBJS)
	    $(CXX) -o $@ $^ $(LIB) $(CXX_FLAGS)

gtest: src/leviosam_test.o src/chain.o src/leviosam_utils.o
	    $(CXX) -o $@ $^ $(LIB) $(CXX_FLAGS) -lgtest_main -lgtest
#gtest: leviosam_test.o chain.o
	    # $(CXX) -o $@ $< $(LIB) $(CXX_FLAGS) -lgtest_main -lgtest

clean:
	    rm *.gch src/*.o $(MAIN) gtest
