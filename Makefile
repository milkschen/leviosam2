CXX=g++
CXX_FLAGS=--std=c++11
INC=
LIB=-lsdsl -lhts

all: liftover

liftover: test.cpp liftover.hpp
	$(CXX) -o $@ $< $(LIB)
