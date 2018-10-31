CXX = g++
CC = g++

CXXFLAGS = -O3 -flto -Wall

all: main

basis_functions.o: basis_functions.cpp basis_functions.h

element.o: element.cpp element.h

linear_element.o: linear_element.cpp linear_element.h

tridiagonal_matrix.o: tridiagonal_matrix.cpp tridiagonal_matrix.h
	$(CXX) $(CXXFLAGS) -framework Accelerate -flax-vector-conversions -c -o $@ $<

utilities.o: utilities.cpp utilities.h

Burgers_model.o: Burgers_model.cpp Burgers_model.h
	$(CXX) $(CXXFLAGS) -fopenmp -c -o $@ $<

main: main.o basis_functions.o element.o linear_element.o tridiagonal_matrix.o utilities.o Burgers_model.o
	$(CXX) $(CXXFLAGS) -framework Accelerate -flax-vector-conversions -fopenmp -o $@ $^

clean:
	rm -f *.o *.dat main
