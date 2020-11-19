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
