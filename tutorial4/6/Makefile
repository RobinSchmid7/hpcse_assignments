.SUFFIXES:
.SUFFIXES: .cpp .o
CXXFLAGS = -O2 -g
CXX = g++
PREFIX = $(HOME)/prefix
O = read.o write.o
L = libmatrix_io.a
H = matrix_io.h
all: $L
$L: $O
	ar r $L $O
.cpp.o:
	$(CXX) $< -c -fPIE $(CXXFLAGS)
install: $L
	mkdir -p $(PREFIX)/include $(PREFIX)/lib && \
	cp $H $(PREFIX)/include/ && \
	cp $L $(PREFIX)/lib/
clean:; -rm $L $O
