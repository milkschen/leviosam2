PRGNAME=leviosam
CXX=g++
CXX_FLAGS=--std=c++11 -lpthread
INC=
LIB=-lsdsl -lhts

all: $(PRGNAME)

leviosam: leviosam.cpp leviosam.hpp
	$(CXX) -o $@ $< $(LIB) $(CXX_FLAGS)

gtest: leviosam_test.cpp leviosam.hpp
	$(CXX) -o $@ $< $(LIB) $(CXX_FLAGS) -lgtest_main -lgtest

clean:
	rm *.o $(PRGNAME)
