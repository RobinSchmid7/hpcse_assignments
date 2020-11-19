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
