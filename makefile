CXX = g++
CC = g++

CXXFLAGS = -O3 -flto -Wall -flax-vector-conversions -fopenmp

LDFLAGS = -framework Accelerate

# cppsrc = $(wildcard src/*.cpp) \
#          $(wildcard src/elements/*.cpp) \
#          $(wildcard src/linear_algebra/*.cpp) \
#          $(wildcard src/utilities/*.cpp)

cppsrc = src/Burgers_model.cpp \
         src/elements/basis_functions.cpp src/elements/element.cpp src/elements/linear_element.cpp \
         src/linear_algebra/linear_algebra.cpp src/linear_algebra/tridiagonal_matrix.cpp \
         src/utilities/utilities.cpp

obj = $(cppsrc:.cpp=.o)

main: src/main.o $(obj)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

test_dirichlet: src/test/test_dirichlet.o $(obj)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

test_periodic: src/test/test_periodic.o $(obj)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

dns: src/thesis/dns.o $(obj)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

clean:
	rm -f main src/main.o \
	      test_dirichlet src/test_dirichlet.o \
	      test_periodic src/test_periodic.o \
	      dns src/thesis/dns.o \
	      $(obj) \
	      data/*.dat \
	      test_data/*.dat \
	      dns_data/*.dat
