CXX=g++
CXX_FLAGS=--std=c++11
INC=
LIB=-lsdsl -lhts

all: liftover

liftover: liftover.cpp liftover.hpp
	$(CXX) -o $@ $< $(LIB)

test: test.cpp liftover.hpp
	$(CXX) -o $@ $< $(LIB)
