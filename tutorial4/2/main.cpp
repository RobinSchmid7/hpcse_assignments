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
