#include <mkl.h>
#include <stdio.h>
int main()
{
    const int N = 2 << 24, iterations = 64, alignment = 64, inc = 1;
    int i, iter;
    double start, manualResult, blasResult, byteCount, intensity, blasTime,
        flopCount, manualTime, *a, *b;
    a = (double*)mkl_malloc(N * sizeof *a, alignment);
    b = (double*)mkl_malloc(N * sizeof *b, alignment);
    if (a == NULL || b == NULL) {
        fprintf(stderr, "p1: mkl_malloc failed\n");
        return 2;
    }
    for (i = 0; i < N; i++) {
        a[i] = (double)i / N;
        b[i] = (double)i / N;
    }
    dsecnd();
    start = dsecnd();
    for (iter = 0; iter < iterations; iter++) {
        manualResult = 0.0;
        for (i = 0; i < N; i++)
            manualResult += a[i] * b[i];
    }
    manualTime = dsecnd() - start;
    start = dsecnd();
    for (iter = 0; iter < iterations; iter++)
        blasResult = cblas_ddot(N, a, inc, b, inc);
    blasTime = dsecnd() - start;
    flopCount = (double)N * iterations * 2; // One multiplication and one sum
    byteCount = (double)N * iterations * 4 *
                sizeof(double); // Three reads and one write
    intensity = flopCount / byteCount;
    printf("Size: %d\n", N);
    printf("MFlops Executed: %.2f\n", flopCount / (1024.0 * 1024.0));
    printf("MBytes Accessed: %.2f\n", byteCount / (1024.0 * 1024.0));
    printf("Operational Intensity: %f\n", intensity);
    printf("Manual:\n");
    printf("    Result: %.16e\n", manualResult);
    printf("    Computation Time: %.3fs\n", manualTime);
    printf("    GFlops/s:  %f\n",
           flopCount / (manualTime * 1024 * 1024 * 1024));
    printf("BLAS:\n");
    printf("    Result: %.16e\n", blasResult);
    printf("    Computation Time: %.3fs\n", blasTime);
    printf("    GFlops/s:  %f\n", flopCount / (blasTime * 1024 * 1024 * 1024));
    printf("Performance Ratio: %.2f%%\n", (blasTime / manualTime) * 100);
    mkl_free(a);
    mkl_free(b);
}
