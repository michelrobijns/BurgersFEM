CXX = g++
CC = g++

CXXFLAGS = -O3 -flto -Wall -flax-vector-conversions -fopenmp
#CXXFLAGS = -g -flax-vector-conversions -fopenmp

LDFLAGS = -framework Accelerate

# cppsrc = $(wildcard src/*.cpp) \
#          $(wildcard src/elements/*.cpp) \
#          $(wildcard src/linear_algebra/*.cpp) \
#          $(wildcard src/utilities/*.cpp)

cppsrc = src/Burgers_model.cpp \
         src/elements/basis_functions.cpp src/elements/element.cpp src/elements/linear_element.cpp \
         src/linear_algebra/linear_algebra.cpp src/linear_algebra/tridiagonal_matrix.cpp \
         src/utilities/utilities.cpp src/utilities/legendre_rule.cpp

obj = $(cppsrc:.cpp=.o)

main: main.o $(obj)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

clean:
	rm -f main src/main.o \
	      test/test_dirichlet test/test_dirichlet.o \
	      test/test_periodic test/test_periodic.o \
	      thesis/dns thesis/dns.o thesis/projector_model.o \
	      $(obj) \
	      data/*.dat \
	      test/test_data/*.dat \
	      thesis/dns_data/*.dat
