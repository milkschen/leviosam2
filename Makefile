PRGNAME=liftover
CXX=g++
CXX_FLAGS=--std=c++11
INC=
LIB=-lsdsl -lhts

all: $(PRGNAME)

liftover: liftover.cpp liftover.hpp
	$(CXX) -o $@ $< $(LIB) $(CXX_FLAGS)

test: test.cpp liftover.hpp
	$(CXX) -o $@ $< $(LIB) $(CXX_FLAGS)

clean:
	rm *.o $(PRGNAME)
