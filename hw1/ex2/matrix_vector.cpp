/*Author: Robin Schmid, 16-927-725, schmirob@ethz.ch*/
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <chrono>
#include <math.h>

std::vector<double> Ax_row( std::vector<double> &A, std::vector<double> &x ){

  size_t N = x.size();
  std::vector<double> y(N);

  // TODO: Question 2a: Matrix vector multiplication with a row major matrix
    for (size_t i=0; i<N; ++i)
        for (size_t j=0; j<N; ++j)
            y[i] += A[i*N+j]*x[j];

  return y;
}



std::vector<double> Ax_col( std::vector<double> &A, std::vector<double> &x ){

  size_t N = x.size();
  std::vector<double> y(N);

  // TODO: Question 2a: Matrix vector multiplication with a column major matrix
    for (size_t i=0; i<N; ++i)
        for (size_t j=0; j<N; ++j)
            y[i] += A[j*N+i]*x[j];

  return y;
}


double benchmark_Ax( std::vector<double> &A, std::vector<double> &x, bool row_major, double Ns){

  double times = 0;

  for( size_t i=0; i<Ns; i++){
    auto t1 = std::chrono::system_clock::now();
    // TODO: Question 2a: Call the function to be benchmarked
    if(row_major){
    Ax_row(A, x);
    }
    else{
    Ax_col(A, x);
    }
    auto t2 = std::chrono::system_clock::now();
    times += std::chrono::duration<double>(t2-t1).count();
  }
  printf("Done in total %9.4fs  --  average %9.4fs\n", times, times/Ns);

  return times/Ns;
}


int main( int argc, char **argv ) {

    if (argc < 3) {
        printf("Usage: %s [N|matrix dimension] [Ns|number of iterations]\n", argv[0]);
        exit(0);
    }

    size_t N = atoi(argv[1]);
    size_t Ns = atoi(argv[2]);
    std::vector<double> A(N * N), B(N * N), x(N);

    // TODO: Question 2a: Initialize matrices and vector
    //       store A as row major and B as column major
    // Store A row major
    for (size_t i=0; i<N; ++i)
        for (size_t j=0; j<N; ++j)
            A[i*N+j] = i+j;

    // Store B column major
    for (size_t i=0; i<N; ++i)
        for (size_t j=0; j<N; ++j)
            B[j*N+i] = i+j;

    // Store x vector
    for (size_t i=0; i<N; ++i)
        x[i] = i;

    printf("Working with matrix of dimension %zu\n", N);

    printf("A*x (row major).\n");
    double times1 = benchmark_Ax(A, x, true, Ns);

    printf("A*x (column major).\n");
    double times2 = benchmark_Ax(B, x, false, Ns);

    printf("-----------------\n");
    printf("Speedup %.8fs\n", times1 / times2);

    return 0;
}
