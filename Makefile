MAIN=leviosam
CXX=g++
CXX_FLAGS=--std=c++11 -lpthread
LIB=-lsdsl -lhts

DEPS = leviosam.hpp chain.hpp
OBJS = leviosam.o chain.o

all: $(MAIN)

%.o: %.cpp $(DEPS)
	$(CXX) -c -o $@ $< $(LIB) $(CXX_FLAGS)

$(MAIN): $(OBJS)
	$(CXX) -o $@ $^ $(LIB) $(CXX_FLAGS)

gtest: leviosam_test.cpp leviosam.hpp
	$(CXX) -o $@ $< $(LIB) $(CXX_FLAGS) -lgtest_main -lgtest

clean:
	rm *.o $(MAIN)
