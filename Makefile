PRGNAME=leviosam
CXX=g++
CXX_FLAGS=--std=c++11 -lpthread
INC=
LIB=-lsdsl -lhts

all: $(PRGNAME)

leviosam: leviosam.cpp leviosam.hpp
	$(CXX) -o $@ $< $(LIB) $(CXX_FLAGS)

test: test.cpp leviosam.hpp
	$(CXX) -o $@ $< $(LIB) $(CXX_FLAGS)

clean:
	rm *.o $(PRGNAME)
