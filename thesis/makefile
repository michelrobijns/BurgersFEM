CXX = g++
CC = g++

CXXFLAGS = -O3 -flto -Wall -flax-vector-conversions -fopenmp

LDFLAGS = -framework Accelerate

cppsrc = ../src/Burgers_model.cpp \
         ../src/elements/basis_functions.cpp ../src/elements/element.cpp ../src/elements/linear_element.cpp \
         ../src/linear_algebra/linear_algebra.cpp ../src/linear_algebra/tridiagonal_matrix.cpp \
         ../src/utilities/utilities.cpp ../src/utilities/legendre_rule.cpp

obj = $(cppsrc:.cpp=.o)

dns: dns.o projector_model.o $(obj)
	$(CXX) $(CXXFLAGS) -o $@ $^ $(LDFLAGS)

clean:
	rm -f dns dns.o projector_model.o \
	      $(obj) \
	      dns_data/*.dat
