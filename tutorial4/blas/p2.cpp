#include <mkl.h>
#include <stdio.h>
double sum(int n, double* a)
{
    int i;
    double ans = 0;
    for (i = 0; i < n; i++)
        ans += a[i];
    return ans;
}

int main()
{
    const int N = 2 << 13, alignment = 64, inc = 1;
    int i, j;
    double start, byteCount, intensity, blasTime, flopCount, manualTime, *A,
        *b, *manualResult, *blasResult;
    A = (double*)mkl_malloc(N * N * sizeof *A, alignment);
    b = (double*)mkl_malloc(N * sizeof *b, alignment);
    manualResult = (double*)mkl_malloc(N * sizeof *manualResult, alignment);
    blasResult = (double*)mkl_malloc(N * sizeof *blasResult, alignment);
    if (A == NULL || b == NULL || manualResult == NULL || blasResult == NULL) {
        fprintf(stderr, "p1: mkl_malloc failed\n");
        return 2;
    }
    for (i = 0; i < N; i++) {
        b[i] = (double)i / N;
        manualResult[i] = 0;
        blasResult[i] = 0;
        for (j = 0; j < N; j++)
            A[i * N + j] = (double)(i + j) / N;
    }
    dsecnd();
    start = dsecnd();
    for (j = 0; j < N; j++)
      for (i = 0; i < N; i++)
            manualResult[i] += b[j] * A[i * N + j];
    manualTime = dsecnd() - start;
    start = dsecnd();
    cblas_dgemv(CblasRowMajor, CblasNoTrans, N, N, 1.0, A, N, b, inc, 0.0,
                blasResult, inc);
    blasTime = dsecnd() - start;
    flopCount = (double)N * N * 2; // One multiplication and one sum
    byteCount = (double)N * N * 4 * sizeof(double); // Three reads and one write
    intensity = flopCount / byteCount;
    printf("Size: %d\n", N);
    printf("MFlops Executed: %.2f\n", flopCount / (1024.0 * 1024.0));
    printf("MBytes Accessed: %.2f\n", byteCount / (1024.0 * 1024.0));
    printf("Operational Intensity: %f\n", intensity);
    printf("Manual:\n");
    printf("    Result: %.16e\n", sum(N, manualResult));
    printf("    Computation Time: %.3fs\n", manualTime);
    printf("    GFlops/s:  %f\n",
           flopCount / (manualTime * 1024 * 1024 * 1024));
    printf("BLAS:\n");
    printf("    Result: %.16e\n", sum(N, blasResult));
    printf("    Computation Time: %.3fs\n", blasTime);
    printf("    GFlops/s:  %f\n", flopCount / (blasTime * 1024 * 1024 * 1024));
    printf("Performance Ratio: %.2f%%\n", (blasTime / manualTime) * 100);
    mkl_free(A);
    mkl_free(b);
    mkl_free(manualResult);
    mkl_free(blasResult);
}
