libmatrix_io.a: lib.o
	ar r libmatrix_io.a lib.o
lib.o: lib.cpp matrix_io.h
	g++ -c -fPIE lib.cpp
install: libmatrix_io.a
	mkdir -p $(HOME)/prefix/include $(HOME)/prefix/lib && \
	cp matrix_io.h $(HOME)/prefix/include/ && \
	cp libmatrix_io.a $(HOME)/prefix/lib/
