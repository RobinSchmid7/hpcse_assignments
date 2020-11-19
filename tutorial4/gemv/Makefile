.SUFFIXES:
.SUFFIXES: .cpp
PREFIX = $(HOME)/prefix
CXXFLAGS = -O2 -g
CXX = g++
MKL_FLAGS = `pkg-config --cflags --libs mkl-static-lp64-seq`
M = gemv
all: $M
gemv: timer.h
.cpp:
	$(CXX) $< -std=c++11 -fopenmp -mavx2 -mfma $(CXXFLAGS) $(LDFLAGS) $(MKL_FLAGS) -o $@
clean:
	rm -f -- $M
