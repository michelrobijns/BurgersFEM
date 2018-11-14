CXX = g++
CC = g++

CXXFLAGS = -O3 -flto -Wall -flax-vector-conversions -fopenmp

LDFLAGS = -framework Accelerate

cppsrc = ../src/Burgers_model.cpp \
         ../src/elements/basis_functions.cpp ../src/elements/element.cpp ../src/elements/linear_element.cpp \
         ../src/linear_algebra/linear_algebra.cpp ../src/linear_algebra/tridiagonal_matrix.cpp \
         ../src/utilities/utilities.cpp ../src/utilities/legendre_rule.cpp

obj = $(cppsrc:.cpp=.o)

all: test_dirichlet test_periodic

test_dirichlet: test_dirichlet.o $(obj)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

test_periodic: test_periodic.o $(obj)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

clean:
	rm -f main src/main.o \
	      test_dirichlet test_dirichlet.o \
	      test_periodic test_periodic.o \
	      $(obj) \
	      test_data/*.dat \
