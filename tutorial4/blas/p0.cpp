#include <cstring>
#include <matrix_io.h>
#include <mkl.h>
#include <stdio.h>
#include <stdlib.h>
#include <vector>

int main(int argc, char** argv)
{
    enum { inc = 1 };
    int ma, na, mb, nb;
    double *a, *b, alpha, beta;
    char* function;
    FILE* file;
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
        if (na != mb || nb != 1) {
            fprintf(stderr, "p0: wrong dimensions for dgemv: %dx%d and %dx%d\n",
                    ma, na, mb, nb);
            return 2;
        }
        std::vector<double> y(na * nb);
        alpha = 1.0;
        beta = 0.0;
        cblas_dgemv(CblasRowMajor, CblasNoTrans, ma, na, alpha, a, na, b, inc,
                    beta, y.data(), inc);
        if (matrix_io_write(stdout, ma, 1, y.data()) != 0)
            return 2;
    } else {
        fprintf(stderr, "p0: unknown function '%s'\n", function);
        return 2;
    }
    free(b);
    free(a);
}
