.SUFFIXES:
.SUFFIXES: .cpp
PREFIX = $(HOME)/prefix
CXXFLAGS = -O2 -g
CXX = g++
MKL_FLAGS = `pkg-config --cflags --libs mkl-static-lp64-seq`
MATRIX_IO_FLAGS = -I$(PREFIX)/include -L$(PREFIX)/lib -lmatrix_io
M = p0 p1 p2
all: $M
.cpp:
	$(CXX) $< $(CXXFLAGS) $(LDFLAGS) $(MKL_FLAGS) $(MATRIX_IO_FLAGS) -o $@
clean:
	rm -f -- $M
