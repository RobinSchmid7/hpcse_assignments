/*Author: Robin Schmid, 16-927-725, schmirob@ethz.ch*/
#include <omp.h>
#include <random>
#include <stdio.h>
#include <stdlib.h>

// Integrand
inline double F(double x, double y)
{
    if (x * x + y * y < 1.) { // inside unit circle
        return 4.;
    }
    return 0.;
}

// Method 0: serial
double C0(size_t n)
{
    // random generator with seed 0
    std::default_random_engine g(0);
    // uniform distribution in [0, 1]
    std::uniform_real_distribution<double> u;

    double s = 0.; // sum
    for (size_t i = 0; i < n; ++i) {
        double x = u(g);
        double y = u(g);
        s += F(x, y);
    }
    printf("Serial: %f\n", s/n);
    return s / n;
}

// Method 1: openmp, no arrays
// TODO: Question 1a.1
double C1(size_t n)
{
    double s = 0.;
    #pragma omp parallel // Spawn a team of threads
    {
        const int tid = omp_get_thread_num();
        std::default_random_engine g(tid); // Initialize each thread with a different random seed
        #pragma omp for reduction(+: s) // Optimize the sum with the reduction directive
        for (size_t i = 0; i < n; ++i) {
            std::uniform_real_distribution<double> u; // Use in for loop for better efficiency
            double x = u(g);
            double y = u(g);
            s += F(x, y);
        }
    }
    printf("C1: %f\n", s/n);
    return s / n;
}

// Method 2, only `omp parallel for reduction`, arrays without padding
// TODO: Question 1a.2
double C2(size_t n)
{
    int nthreads;
    nthreads = omp_get_max_threads();
    double s = 0.;
    // Need to generate the random numbers outside of the for loop otherwise not random
    // Use vector with individual random generators for each thread
    std::vector<std::default_random_engine> g(nthreads);
    for (size_t i = 0; i < nthreads; ++i) {
        g[i].seed(2*i+1); // Individual seed for each thread, use linear scaling s.t. seeds are different and not 0
    }

    #pragma omp parallel for reduction(+: s) // Only allowable directive
    for (size_t i = 0; i < n; ++i) {
        std::uniform_real_distribution<double> u;
        const int tid = omp_get_thread_num();
        double x = u(g[tid]); // Use individual generated seed for each thread
        double y = u(g[tid]);
        s += F(x, y);

    }
    return s / n;
}

// Method 3, only `omp parallel for reduction`, arrays with padding
// TODO: Question 1a.3
double C3(size_t n)
{
    int nthreads;
    nthreads = omp_get_max_threads();
    double s = 0.;
    /* Need to generate the random numbers outside of the for loop otherwise not random
     Each generator works on its own cache line of size of 64 byte on Euler to prevent race conditions
     Maximal sizeof(std::default_random_engine): 2496 byte
     Actual sizeof(std::default_random_engine) varies depending on the value, so if using
     a fixed array size a lot of memory is wasted since we use the maximal size of this random type*/

     /* Option 1: Use spacing with maximal size of array, works but too much memory for spacing
     since some sizeof(std::default_random_engine) << 2496 byte*/
     //std::vector<std::default_random_engine> g(64 * nthreads));
     /*Option 2: Use a struct to space individual random generators*/
    struct g_struct {
        std::default_random_engine g;
        double space[8]; // Use spacing of 64 byte (8*8 double), ensures a new cache line
    };

    std::vector<g_struct> g_pad(nthreads); // Padded array

    for (size_t i = 0; i < nthreads; ++i) {
        g_pad[i].g.seed(2*i+1); // Individual seed for each thread, use linear scaling s.t. seeds are different and not 0
    }

#pragma omp parallel for reduction(+: s) // Only allowed
    for (size_t i = 0; i < n; ++i) {
        std::uniform_real_distribution<double> u;
        const int tid = omp_get_thread_num();
        double x = u(g_pad[tid].g); // Use individual generated seed for each thread
        double y = u(g_pad[tid].g);
        s += F(x, y);

    }
    return s / n;
}

// Another option adapted from lecture slides, here use additional directives which are not allowed for this exercise
/*
double C3(size_t n)
{
    // Create an array where each thread can store its result
    int nthreads;
    #pragma omp parallel
    #pragma omp master
    nthreads = omp_get_max_threads();
    // Use an array in which each thread can store its results, use padding s.t.
    // there are now cache lines used by different threads at the same time, in
    // a cache line of size 64 can place 8 doubles
    double* s = new double[8*nthreads];

    #pragma omp parallel
    {
        const int tid = omp_get_thread_num();
        // Initialize each thread with a different seed
        std::default_random_engine g(tid);
        std::uniform_real_distribution<double> u;
        #pragma omp for nowait
        for (size_t i = 0; i < n; ++i) {
            double x = u(g);
            double y = u(g);
            s[8*tid] += F(x, y);
        }
    }
    // Sum over all results of the individual threads
    for (size_t i = 1; i < nthreads; ++i)
        s[0] += s[8*i];

//delete[] s; // removes effect of padding?
    printf("C3: %f\n", s[0]/n);
    return s[0] / n;
}*/

// Returns integral of F(x,y) over unit square (0 < x < 1, 0 < y < 1).
// n: number of samples
// m: method
double C(size_t n, size_t m)
{
    switch (m) {
    case 0:
        return C0(n);
    case 1:
        return C1(n);
    case 2:
        return C2(n);
    case 3:
        return C3(n);
    default:
        printf("Unknown method '%ld'\n", m);
        abort();
    }
}

int main(int argc, char* argv[])
{
    // default number of samples
    const size_t ndef = 1e8;

    if (argc < 2 || argc > 3 || std::string(argv[1]) == "-h") {
        fprintf(stderr, "usage: %s METHOD [N=%ld]\n", argv[0], ndef);
        fprintf(stderr, "Monte-Carlo integration with N samples.\n\
METHOD:\n\
0: serial\n\
1: openmp, no arrays\n\
2: `omp parallel for reduction`, arrays without padding\n\
3: `omp parallel for reduction`, arrays with padding\n");
        return 1;
    }

    // method
    size_t m = atoi(argv[1]);
    // number of samples
    size_t n = (argc > 2 ? atoi(argv[2]) : ndef);
    // reference solution
    double ref = 3.14159265358979323846;

    double wt0 = omp_get_wtime();
    double res = C(n, m);
    double wt1 = omp_get_wtime();

    printf("res:  %.20f\nref:  %.20f\nerror: %.20e\ntime: %.20f\n", res, ref,
           res - ref, wt1 - wt0);

    return 0;
}
