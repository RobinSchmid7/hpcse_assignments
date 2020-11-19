#include <functional>
#include <immintrin.h>
#include <mkl.h>
#include <stdlib.h>
#include <string.h>
#include <string>

#include "timer.h"

// AVX
const size_t simd_width_bytes = 256 / 8;

static double flush_cache()
{
    const size_t n = 10 * 1024 * 1024;
    static double* junk = nullptr;
    if (!junk)
        junk = new double[n];

    volatile double tmp = 0;
    for (size_t i = 0; i < n; ++i)
        junk[i] = i;
    for (size_t i = 0; i < n; ++i)
        tmp += junk[i];
    return tmp;
}

template <typename real>
void initialise(int m, int n, int lda, real* A, real* x)
{
    memset(A, 0, lda * n * sizeof(real));
    for (int i = 0; i < m; ++i) {
        x[i] = 2;
        for (int j = 0; j < n; ++j)
            A[i * lda + j] = 1;
    }
}

template <typename real> real check(int n, const real* y)
{
    real sum = 0.f;
    for (int i = 0; i < n; ++i)
        sum += y[i];
    volatile real res = sum / n;
    return res;
}

template <typename real>
void gemv_naive(int m, int n, int lda, const real* A, const real* x, real* y)
{
    for (int i = 0; i < m; ++i)
        y[i] = 0;

    for (int i = 0; i < m; ++i)
        for (int j = 0; j < n; ++j)
            y[i] += A[i * lda + j] * x[j];
}

template <typename real>
void gemv_block2(int m, int n, int lda, const real* A, const real* x, real* y)
{
    real y0, y1, xj;

    for (int i = 0; i < m; i += 2) {
        y0 = y1 = 0;
        for (int j = 0; j < n; ++j) {
            xj = x[j];
            y0 += A[(i + 0) * lda + j] * xj;
            y1 += A[(i + 1) * lda + j] * xj;
        }
        y[i + 0] = y0;
        y[i + 1] = y1;
    }
}

template <typename real>
void gemv_block4(int m, int n, int lda, const real* A, const real* x, real* y)
{
    real y0, y1, y2, y3, xj, a0, a1, a2, a3;

    for (int i = 0; i < m; i += 4) {
        y0 = y1 = y2 = y3 = 0;
        for (int j = 0; j < n; ++j) {
            xj = x[j];
            a0 = A[(i + 0) * lda + j];
            a1 = A[(i + 1) * lda + j];
            a2 = A[(i + 2) * lda + j];
            a3 = A[(i + 3) * lda + j];

            y0 += a0 * xj;
            y1 += a1 * xj;
            y2 += a2 * xj;
            y3 += a3 * xj;
        }
        y[i + 0] = y0;
        y[i + 1] = y1;
        y[i + 2] = y2;
        y[i + 3] = y3;
    }
}

template <typename real>
void gemv_block_omp(int m, int n, int lda, const real* A, const real* x,
                    real* y)
{
    real y0, y1, xj;

#pragma omp parallel for private(y0, y1, xj)
    for (int i = 0; i < m; i += 2) {
        y0 = y1 = 0;
        for (int j = 0; j < n; ++j) {
            xj = x[j];
            y0 += A[(i + 0) * lda + j] * xj;
            y1 += A[(i + 1) * lda + j] * xj;
        }
        y[i + 0] = y0;
        y[i + 1] = y1;
    }
}

void gemv_simd_s(int m, int n, int lda, const float* A, const float* x,
                 float* y)
{
    int i, j;
    __m256 AA, xx, yy;

    memset(y, 0, m * sizeof(float));
    const int simd_width = simd_width_bytes / sizeof(float);

    for (i = 0; i < m; ++i) {
        yy = _mm256_set_ps(0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f);

        for (j = 0; j < n; j += simd_width) {
            AA = _mm256_load_ps(A + i * lda + j);
            xx = _mm256_load_ps(x + j);
            yy = _mm256_fmadd_ps(AA, xx, yy);
        }

        y[i] += ((yy[0] + yy[1]) + (yy[2] + yy[3])) +
                ((yy[4] + yy[5]) + (yy[6] + yy[7]));
    }
}

void gemv_simd_s_block2(int m, int n, int lda, const float* A, const float* x,
                        float* y)
{
    int i, j;
    __m256 AA0, AA1, xx, yy0, yy1;

    memset(y, 0, m * sizeof(float));
    const int simd_width = simd_width_bytes / sizeof(float);

    for (i = 0; i < m; i += 2) {
        yy0 = yy1 = _mm256_set_ps(0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f);

        for (j = 0; j < n; j += simd_width) {
            AA0 = _mm256_load_ps(A + (i + 0) * lda + j);
            AA1 = _mm256_load_ps(A + (i + 1) * lda + j);
            xx = _mm256_load_ps(x + j);

            yy0 = _mm256_fmadd_ps(AA0, xx, yy0);
            yy1 = _mm256_fmadd_ps(AA1, xx, yy1);
        }

        y[i + 0] += ((yy0[0] + yy0[1]) + (yy0[2] + yy0[3])) +
                    ((yy0[4] + yy0[5]) + (yy0[6] + yy0[7]));
        y[i + 1] += ((yy1[0] + yy1[1]) + (yy1[2] + yy1[3])) +
                    ((yy1[4] + yy1[5]) + (yy1[6] + yy1[7]));
    }
}

void gemv_simd_s_block4(int m, int n, int lda, const float* A, const float* x,
                        float* y)
{
    int i, j;
    __m256 AA0, AA1, AA2, AA3, xx, yy0, yy1, yy2, yy3;

    memset(y, 0, m * sizeof(float));
    const int simd_width = simd_width_bytes / sizeof(float);

    for (i = 0; i < m; i += 4) {
        yy0 = yy1 = yy2 = yy3 =
            _mm256_set_ps(0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f);

        for (j = 0; j < n; j += simd_width) {
            AA0 = _mm256_load_ps(A + (i + 0) * lda + j);
            AA1 = _mm256_load_ps(A + (i + 1) * lda + j);
            AA2 = _mm256_load_ps(A + (i + 2) * lda + j);
            AA3 = _mm256_load_ps(A + (i + 3) * lda + j);
            xx = _mm256_load_ps(x + j);

            yy0 = _mm256_fmadd_ps(AA0, xx, yy0);
            yy1 = _mm256_fmadd_ps(AA1, xx, yy1);
            yy2 = _mm256_fmadd_ps(AA2, xx, yy2);
            yy3 = _mm256_fmadd_ps(AA3, xx, yy3);
        }

#define REDUCE_8(yy)                                                           \
    ((yy[0] + yy[1]) + (yy[2] + yy[3])) + ((yy[4] + yy[5]) + (yy[6] + yy[7]));

        y[i + 0] += REDUCE_8(yy0);
        y[i + 1] += REDUCE_8(yy1);
        y[i + 2] += REDUCE_8(yy2);
        y[i + 3] += REDUCE_8(yy3);

#undef REDUCE_8
    }
}

void gemv_simd_d(int m, int n, int lda, const double* A, const double* x,
                 double* y)
{
    int i, j;
    __m256d AA, xx, yy;

    memset(y, 0, m * sizeof(double));
    const int simd_width = simd_width_bytes / sizeof(double);

    for (i = 0; i < m; ++i) {
        yy = _mm256_set_pd(0., 0., 0., 0.);

        for (j = 0; j < n; j += simd_width) {
            AA = _mm256_load_pd(A + i * lda + j);
            xx = _mm256_load_pd(x + j);
            yy = _mm256_fmadd_pd(AA, xx, yy);
        }

        y[i] += (yy[0] + yy[1]) + (yy[2] + yy[3]);
    }
}

void gemv_simd_d_block2(int m, int n, int lda, const double* A, const double* x,
                        double* y)
{
    int i, j;
    __m256d AA0, AA1, xx, yy0, yy1;

    memset(y, 0, m * sizeof(double));
    const int simd_width = simd_width_bytes / sizeof(double);

    for (i = 0; i < m; i += 2) {
        yy0 = yy1 = _mm256_set_pd(0., 0., 0., 0.);

        for (j = 0; j < n; j += simd_width) {
            AA0 = _mm256_load_pd(A + (i + 0) * lda + j);
            AA1 = _mm256_load_pd(A + (i + 1) * lda + j);
            xx = _mm256_load_pd(x + j);

            yy0 = _mm256_fmadd_pd(AA0, xx, yy0);
            yy1 = _mm256_fmadd_pd(AA1, xx, yy1);
        }

        y[i + 0] += (yy0[0] + yy0[1]) + (yy0[2] + yy0[3]);
        y[i + 1] += (yy1[0] + yy1[1]) + (yy1[2] + yy1[3]);
    }
}

void sgemv_blas(int m, int n, int lda, const float* A, const float* x, float* y)
{
    cblas_sgemv(CblasRowMajor, CblasNoTrans, m, n, 1.f, A, lda, x, 1.f, 0.f, y,
                1);
}

void dgemv_blas(int m, int n, int lda, const double* A, const double* x,
                double* y)
{
    cblas_dgemv(CblasRowMajor, CblasNoTrans, m, n, 1.0, A, lda, x, 1.0, 0.0, y,
                1);
}

template <typename real>
using GemvFunc =
    std::function<void(int, int, int, const real*, const real*, real*)>;

template <typename real>
static void time_function(GemvFunc<real> gemv, std::string name)
{
    FILE* f = fopen((name + ".dat").c_str(), "w");

    fprintf(stderr, "profiling %s...\n", name.c_str());

    for (int i = 100; i < 5000; i += 100) {
        int m = i, n = i, lda;
        real *A, *x, *y;
        int simd_width = simd_width_bytes / sizeof(real);
        lda = (m + simd_width - 1) / simd_width;
        lda *= simd_width;

        posix_memalign((void**)&A, simd_width_bytes, lda * n * sizeof(real));
        posix_memalign((void**)&x, simd_width_bytes, m * sizeof(real));
        posix_memalign((void**)&y, simd_width_bytes, n * sizeof(real));

        initialise(m, n, lda, A, x);

        const int nsteps = i < 1000 ? 100 : 20;

        mTimer timer;
        double t = 0;
        for (int step = 0; step < nsteps; ++step) {
            flush_cache();
            timer.start();
            gemv(m, n, lda, A, x, y);
            t += timer.elapsed();
            check(n, y);
        }

        double time = t / nsteps;
        long FLOP = n * m * 2;
        double FLOPS = FLOP / (time * 1e-3);

        fprintf(f, "%d\t%g\t%g\n", i, time, FLOPS);

        free(A);
        free(x);
        free(y);
    }
    fclose(f);
}

int main()
{
    time_function<float>(&gemv_naive<float>, "naive_single");
    time_function<float>(&gemv_block2<float>, "block2_single");
    time_function<float>(&gemv_block4<float>, "block4_single");
    time_function<float>(&gemv_block_omp<float>, "block2+omp_single");
    time_function<float>(&gemv_simd_s, "simd_single");
    time_function<float>(&gemv_simd_s_block2, "simd+block2_single");
    time_function<float>(&sgemv_blas, "blas_single");

    time_function<double>(&gemv_naive<double>, "naive_double");
    time_function<double>(&gemv_block2<double>, "block2_double");
    time_function<double>(&gemv_block4<double>, "block4_double");
    time_function<double>(&gemv_block_omp<double>, "block2+omp_double");
    time_function<double>(&gemv_simd_d, "simd_double");
    time_function<double>(&gemv_simd_d_block2, "simd+block2_double");
    time_function<double>(&dgemv_blas, "blas_double");

    return 0;
}
