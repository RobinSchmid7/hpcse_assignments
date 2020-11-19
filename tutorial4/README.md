I want to make a library `matrix_io` with a function `matrix_io_write`. Here is a prototype of what `matrix_io_write` should do:

```c++
$ cat 1/main.cpp
#include <stdio.h>
int main()
{
    int m = 2, n = 3;
    double a[] = {1, -2, 3, -4, 5, -6};
    int i, j;
    for (i = 0; i < m; i++) {
        for (j = 0; j < n; j++)
            printf("%-+.16e ", a[i * n + j]);
        printf("\n");
    }
}
```

```
$ g++ 1/main.cpp
$ ./a.out
+1.0000000000000000e+00 -2.0000000000000000e+00 +3.0000000000000000e+00
-4.0000000000000000e+00 +5.0000000000000000e+00 -6.0000000000000000e+00
```

Even in this simple case I have to make several design decisions: how to represent a matrix, should I print only to stdout, should printf specifier (`%-+.16e `) be a parameter, what to do on error. I settled on this:

```c++
$ cat 2/main.cpp
#include <stdio.h>
int matrix_io_write(FILE* file, int m, int n, double* a)
{
    int i, j;
    for (i = 0; i < m; i++) {
        for (j = 0; j < n; j++) {
            if (fprintf(file, "%-+.16e ", a[i * n + j]) < 0)
                return 1;
        }
        if (fprintf(file, "\n") < 0)
            return 1;
    }
    return 0;
}

int main()
{
    int m = 2, n = 3;
    double a[] = {1, -2, 3, -4, 5, -6};
    matrix_io_write(stdout, m, n, a);
}
```

Here the code of the library `matrix_io_write` and of the client `main` are in the same file. In the following I split them:

```c++
$ cat 3/lib.cpp
#include <stdio.h>
int matrix_io_write(FILE*, int, int, double*);
int matrix_io_write(FILE* file, int m, int n, double* a)
{
    int i, j;
    for (i = 0; i < m; i++) {
        for (j = 0; j < n; j++) {
            if (fprintf(file, "%-+.16e ", a[i * n + j]) < 0)
                return 1;
        }
        if (fprintf(file, "\n") < 0)
            return 1;
    }
    return 0;
}
```

and

```c++
$ cat 3/client.cpp
#include <stdio.h>
int matrix_io_write(FILE*, int, int, double*);
int main()
{
    int m = 2, n = 3;
    double a[] = {1, -2, 3, -4, 5, -6};
    if (matrix_io_write(stdout, m, n, a) != 0) {
        fprintf(stderr, "client: matrix_io_write failed\n");
        return 1;
    }
}
```

Now the library can be compiled separatly and linked with the client.

```
$ g++ -c 3/lib.cpp
$ g++ 3/client.cpp lib.o
$ ./a.out
+1.0000000000000000e+00 -2.0000000000000000e+00 +3.0000000000000000e+00
-4.0000000000000000e+00 +5.0000000000000000e+00 -6.0000000000000000e+00
```

`matrix_io_write` must be declared in the client. I put the declaration in
a header file `matrix_io.h` which the client includes. A library can also
include `matrix_io.h` to be sure the declaration is consistent with
the definition.

```c++
$ cat 4/matrix_io.h
int matrix_io_write(FILE*, int, int, double *);
$ cat 4/lib.cpp
#include "matrix_io.h"
#include <stdio.h>
int matrix_io_write(FILE* file, int m, int n, double* a)
{
    int i, j;
    for (i = 0; i < m; i++) {
        for (j = 0; j < n; j++) {
            if (fprintf(file, "%-+.16e ", a[i * n + j]) < 0)
                return 1;
        }
        if (fprintf(file, "\n") < 0)
            return 1;
    }
    return 0;
}
```

Since there will be several steps with dependencies I create a makefile

```
$ cat 4/Makefile
libmatrix_io.a: lib.o
        ar r libmatrix_io.a lib.o
lib.o: lib.cpp matrix_io.h
        g++ -c lib.cpp
install: libmatrix_io.a
        mkdir -p $(HOME)/prefix/include $(HOME)/prefix/lib && \
        cp matrix_io.h $(HOME)/prefix/include && \
        cp libmatrix_io.a $(HOME)/prefix/lib
```
It builds an objects file `lib.o` and put it into an archive `libmatrix_io.a`. Usually an archive contains several object files. `install` target copies files to a location where it is convenient for a client to find them.
```
$ (cd 4 && make install)
g++ -c lib.cpp
ar r libmatrix_io.a lib.o
ar: creating libmatrix_io.a
mkdir -p /home/lisergey/prefix/include /home/lisergey/prefix/lib && \
cp matrix_io.h /home/lisergey/prefix/include && \
cp libmatrix_io.a /home/lisergey/prefix/lib
```
To compile a client with the installed library the locations of the
header and of the archive  must be set, `-lmatrix_io` is an instruction to
search for a file `libmatrix_io.a`:
```c++
$ cat 5/client.cpp
#include <matrix_io.h>
#include <stdio.h>
int main()
{
    int m = 2, n = 3;
    double a[] = {1, -2, 3, -4, 5, -6};
    if (matrix_io_write(stdout, m, n, a) != 0) {
        fprintf(stderr, "client: matrix_io_write failed\n");
        return 1;
    }
}

$ g++ -I$HOME/prefix/include -L$HOME/prefix/lib 5/client.cpp -lmatrix_io
$ ./a.out
+1.0000000000000000e+00 -2.0000000000000000e+00 +3.0000000000000000e+00
-4.0000000000000000e+00 +5.0000000000000000e+00 -6.0000000000000000e+00
```

Finally, I add a function `matrix_io_read()` to `6/read.cpp` and put
`matrix_io_write()` into `6/write.cpp`, I also create a more flexible
Makefile:

```
$ cat 6/matrix_io.h
#include <stdio.h>
int matrix_io_write(FILE*, int, int, double*);
int matrix_io_read(FILE*, int*, int*, double**);
$ cat 6/Makefile
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
        $(CXX) $< -c $(CXXFLAGS)
install: $L
        mkdir -p $(PREFIX)/include $(PREFIX)/lib && \
        cp $H $(PREFIX)/include && \
        cp $L $(PREFIX)/lib
clean:; -rm $L $O
```
`CXX`, `CXXFLAGS`, 'PREFIX` are parameters. For example,
```
$ (cd 6 && make 'CXX = clang++' 'CXXFLAGS = -Wall -Wextra -O3 -g')
clang++ read.cpp -c -Wall -Wextra -O3 -g
clang++ write.cpp -c -Wall -Wextra -O3 -g
ar r libmatrix_io.a read.o write.o
ar: creating libmatrix_io.a
```
Install a library
```
$ (cd 6 && make clean install)
```
[blas/Makefile](blas/Makefile) compiles programs which use `matrix_io` and Intel MKL. A program `pkg-config` provides Intel MKL compilation flags.
```
$ cat blas/Makefile
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

```

I put several test matrices in files
```
$ cat blas/a
10, 20, 30
$ cat blas/b
1
2
3
$ cat blas/Am
1,2,3
4,5,6
$ cat blas/x
10
20
30
```

A program [blas/p0.cpp](blas/p0.cpp) takes a name of a function (dot or gemv) and two matrices and calls a corresponding MKL function.

```c++
$ cat blas/p0.cpp
#include <stdio.h>
#include <stdlib.h>
#include <cstring>
#include <vector>
#include <mkl.h>
#include <matrix_io.h>

int main(int argc, char** argv)
{
  enum {inc = 1};
  int ma, na, mb, nb;
  double *a, *b, alpha, beta;
  char *function;
  FILE *file;
  if (argc < 4) {
    fprintf(stderr, "p0: needs three arguments\n");
    return 1;
  }
  argv++;
  function = *argv++;
  if ((file = fopen(*argv, "r")) == NULL ||
      matrix_io_read(file, &ma, &na, &a) != 0) {
    fprintf(stderr, "p0: fail to read '%s'\n", *argv);
  }
  argv++;
  fclose(file);
  if ((file = fopen(*argv, "r")) == NULL ||
      matrix_io_read(file, &mb, &nb, &b) != 0) {
    fprintf(stderr, "p0: fail to read '%s'\n", *argv);
  }
  fclose(file);

  if (strcmp(function, "dot") == 0) {
    if (ma != 1 || na != mb || nb != 1) {
      fprintf(stderr, "p0: wrong dimensions for dot: %dx%d and %dx%d\n",
              ma, na, mb, nb);
      return 2;
    }
    printf("%g\n", cblas_ddot(na, a, inc, b, inc));
  } else if (strcmp(function, "gemv") == 0) {
    if (na != mb  || nb != 1) {
      fprintf(stderr, "p0: wrong dimensions for dgemv: %dx%d and %dx%d\n",
              ma, na, mb, nb);
      return 2;
    }
    std::vector<double> y(na * nb);
    alpha = 1.0;
    beta = 0.0;
    cblas_dgemv(CblasRowMajor, CblasNoTrans, ma, na, alpha,
                a, na, b, inc, beta, y.data(), inc);
    if (matrix_io_write(stdout, ma, 1, y.data()) != 0)
      return 2;
  } else {
    fprintf(stderr, "p0: unknown function '%s'\n", function);
    return 2;
  }
  free(b);
  free(a);
  mkl_finalize();
}
```

For example,
```
$ cd blas
$ module load new mkl/2018.1
$ make
g++ p0.cpp -O2 -g  `pkg-config --cflags --libs mkl-static-lp64-seq` -I/cluster/home/lisergey/prefix/include -L/cluster/home/lisergey/prefix/lib -lmatrix_io -o p0
g++ p1.cpp -O2 -g  `pkg-config --cflags --libs mkl-static-lp64-seq` -I/cluster/home/lisergey/prefix/include -L/cluster/home/lisergey/prefix/lib -lmatrix_io -o p1
g++ p2.cpp -O2 -g  `pkg-config --cflags --libs mkl-static-lp64-seq` -I/cluster/home/lisergey/prefix/include -L/cluster/home/lisergey/prefix/lib -lmatrix_io -o p2
$ ./p0 dot a b
140
```

If MKL_VERBOSE is set to 1 the calls of MKL functions are traced:

```
$ MKL_VERBOSE=1 ./p0 gemv Am x
MKL_VERBOSE Intel(R) MKL 2020.0 Product build 20191122 for Intel(R) 64 architecture Intel(R) Streaming SIMD Extensions 4.2 (Intel(R) SSE4.2) enabled processors, Lnx 1.60GHz lp64 sequential
MKL_VERBOSE DGEMV(T,3,2,0x7fffb2f9c980,0x55b2fcafc160,3,0x55b2fcafc0a0,1,0x7fffb2f9c988,0x55b2fcafc0c0,1) 53.50us CNR:OFF Dyn:1 FastMM:1 TID:0  NThr:1
+1.4000000000000000e+02
+3.2000000000000000e+02
```

The programs [blas/p1.cpp](blas/p1.cpp) and [blas/p2.cpp](blas/p2.cpp) compare performance of MKL functions with the manual implementations.

```
$ ./p2
Size: 16384
MFlops Executed: 512.00
MBytes Accessed: 8192.00
Operational Intensity: 0.062500
Manual:
    Result: 1.5657096575000000e+08
    Computation Time: 3.345s
    GFlops/s:  0.149491
BLAS:
    Result: 1.5657096575000000e+08
    Computation Time: 0.157s
    GFlops/s:  3.180541
Performance Ratio: 4.70%
```
One can try different compilers and flags. For example,
```
$ module load llvm
$ make clean
$ make 'CXXFLAGS = -O3 -march=native' 'CXX = clang++'
```
or
```
$ module load intel
$ make clean
$ make 'CXXFLAGS = -O3 -march=native' 'CXX = icpc'
```
also consider Link Line Advisor. A program in [gemv/gemv.cpp](gemv/gemv.cpp) compete with Intel MKL using different approaches.

# References
- [intel:mkl_malloc](https://software.intel.com/content/www/us/en/develop/documentation/mkl-developer-reference-c/top/support-functions/memory-management/mkl-malloc.html)
- [intel:cblas_ddot](https://software.intel.com/content/www/us/en/develop/documentation/mkl-developer-reference-c/top/blas-and-sparse-blas-routines/blas-routines/blas-level-1-routines-and-functions/cblas-dot.html)
- [intel:cblas_dgemv](https://software.intel.com/content/www/us/en/develop/documentation/mkl-developer-reference-c/top/blas-and-sparse-blas-routines/blas-routines/blas-level-2-routines/cblas-gemv.html)
- [intel:cblas_dgemm](https://software.intel.com/content/www/us/en/develop/documentation/mkl-developer-reference-c/top/blas-and-sparse-blas-routines/blas-routines/blas-level-3-routines/cblas-gemm.html)
- [intel:dsecnd](https://software.intel.com/content/www/us/en/develop/documentation/mkl-developer-reference-fortran/top/support-functions/timing/second-dsecnd.html)
- [Link Line Advisor](https://software.intel.com/content/www/us/en/develop/articles/intel-mkl-link-line-advisor.html)
- `$ man pkg-config`
