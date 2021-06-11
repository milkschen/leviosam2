MAIN=leviosam
CXX=g++
CXX_FLAGS=--std=c++11 -lpthread
LIB=-lsdsl -lhts -lz
CFLAGS=-g -Wall -O2
CLIB=-lz -lhts

OBJS = leviosam.o chain.o bam_aux.o bam_md.o
CDEPS= bam.h
DEPS = leviosam.hpp chain.hpp

all: $(MAIN)

.c.o: %.c $(CDEPS)
	    $(CC) -c -o $@ $< $(CLIB) $(CFLAGS)

.cpp.o: %.cpp $(DEPS)
	    $(CXX) -c -o $@ $< $(LIB) $(CXX_FLAGS)

$(MAIN): $(OBJS)
	    $(CXX) -o $@ $^ $(LIB) $(CXX_FLAGS)

gtest: leviosam_test.o chain.o
	    $(CXX) -o $@ $^ $(LIB) $(CXX_FLAGS) -lgtest_main -lgtest
	    # $(CXX) -o $@ $< $(LIB) $(CXX_FLAGS) -lgtest_main -lgtest

clean:
	    rm *.o $(MAIN)
