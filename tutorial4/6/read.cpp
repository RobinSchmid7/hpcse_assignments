#include <cstdio>
#include <cstring>
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include "matrix_io.h"

int matrix_io_read(FILE* file, int* pm, int* pn, double** pa)
{
    char line[999];
    char *nxt, *cur;
    double *a, x;
    int m, n, n0;
    std::vector<double> buf;
    for (m = 0; fgets(line, sizeof line, file) != NULL; m++) {
        for (cur = line, n = 0; *cur != '\n' && *cur != '\0';)
            switch (*cur) {
            case ' ':
            case ',':
            case '\t':
            case '\r':
                cur++;
                break;
            default:
                x = strtod(cur, &nxt);
                if (cur == nxt) {
                    fprintf(stderr, "svd: line %d: not a double\n", m + 1);
                    fprintf(stderr, "svd: '%s'\n", cur);
                    return 1;
                }
                n++;
                cur = nxt;
                buf.push_back(x);
                break;
            }
        if (m == 0)
            n0 = n;
        else {
            if (n != n0) {
                fprintf(stderr, "svd: line %d: expect %d numbers, got %d\n",
                        m + 1, n0, n);
                return 1;
            }
        }
    }
    if ((a = (double*)malloc(n * m * sizeof *a)) == NULL) {
        fprintf(stderr, "svd: alloc failed\n");
        return 2;
    }
    memcpy(a, buf.data(), n * m * sizeof *a);
    *pm = m;
    *pn = n;
    *pa = a;
    return 0;
}
