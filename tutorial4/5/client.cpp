#include <stdio.h>
#include <matrix_io.h>
int main()
{
    int m = 2, n = 3;
    double a[] = {1, -2, 3, -4, 5, -6};
    if (matrix_io_write(stdout, m, n, a) != 0) {
        fprintf(stderr, "client: matrix_io_write failed\n");
        return 1;
    }
}
