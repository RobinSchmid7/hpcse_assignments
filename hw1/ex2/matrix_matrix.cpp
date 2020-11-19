/*Author: Robin Schmid, 16-927-725, schmirob@ethz.ch*/
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <chrono>
#include <math.h>
#include <string>


void AB( std::vector<double> &A, std::vector<double> &B ){

  size_t N = sqrt(A.size());
  std::vector<double> C(N*N);

  // TODO: Question 2c: Straightforward matrix-matrix multiplication
  for (size_t i=0; i<N; ++i)
    for (size_t j=0; j<N; ++j)
      for (size_t k=0; k<N; ++k)
        C[i*N+j] += A[i*N+k]*B[k*N+j];
}


void AB_block_row( std::vector<double> &A, std::vector<double> &B, size_t blockSize ){

  size_t N = sqrt(A.size());
  std::vector<double> C(N*N);
  size_t b = blockSize;

  // TODO: Question 2c: Block matrix-matrix multiplication - B in row major
  // Iterate through a line (row / column) of the blocks of the matrix
  for (size_t k=0; k<N/b; ++k)
      // Iterating through all blocks of the matrix
      for (size_t i=0; i<N/b; ++i)
          for (size_t j=0; j<N/b; ++j)
              // Iterate through all elements within a block
              for (size_t bi=0; bi<b; ++bi)
                  for (size_t bj=0; bj<b; ++bj)
                      // Iterating through a line (row / column) of a block
                      for (size_t bk=0; bk<b; ++bk)
                          C[i*N*b+j*b+bi*N+bj] +=
                                  A[i*N*b+k*b+bi*N+bk] *
                                  B[k*N*b+j*b+bk*N+bj];
}


void AB_block_col( std::vector<double> &A, std::vector<double> &B, size_t blockSize ){

  size_t N = sqrt(A.size());
  std::vector<double> C(N*N);
  size_t b = blockSize;
  // TODO: Question 2c: Block matrix-matrix multiplication - B in column major
    // Iterate through a line (row / column) of the blocks of the matrix
    for (size_t k=0; k<N/b; ++k)
        // Iterating through all blocks of the matrix
        for (size_t i=0; i<N/b; ++i)
            for (size_t j=0; j<N/b; ++j)
                // Iterate through all elements within a block
                for (size_t bi=0; bi<b; ++bi)
                    for (size_t bj=0; bj<b; ++bj)
                        for (size_t bk=0; bk<b; ++bk)
                            // Iterating through a line (row / column) of a block
                            C[i*N*b+j*b+bi*N+bj] +=
                                    A[i*N*b+k*b+bi*N+bk] *
                                    B[j*N*b+k*b+bj*N+bk];
}


double benchmark_AB( std::vector<double> &A, std::vector<double> &B, size_t mode, size_t blockSize, size_t Ns ){

  size_t N = sqrt(A.size());
  double times = 0;

  // TODO: Check that the matrix size is divided by the blockSize when mode==2 or 3
  if( (mode==2 or mode==3) &&  N%blockSize!=0 ){
    printf("Error: the size of the matrix (%zu) should be divided by the blockSize variable (%zu).\n",N,blockSize);
    exit(1);
  }

  for( size_t i=0; i<Ns; i++){
    auto t1 = std::chrono::system_clock::now();
    // TODO: Question 2c: Call the function to be benchmarked
    if( mode==1 ){
      AB(A, B);
    }
    else if( mode==2 ){
      AB_block_row(A, B, blockSize);
    }
    else if( mode==3 ){
      AB_block_col(A, B, blockSize);
    }
    auto t2 = std::chrono::system_clock::now();
    times += std::chrono::duration<double>(t2-t1).count();
  }
  printf("Done in total %9.4fs  --  average %9.4fs\n", times, times/Ns);

  return times/Ns;

}


int main(){

    std::vector<int> matrixSize{256, 512, 1024, 2048};
    size_t M = matrixSize.size();

    std::vector<size_t> blockSize{2, 4, 8, 16, 32, 64, 128};
    size_t Bs = blockSize.size();

    size_t Ns = 5;

    std::vector<double> times1(M);
    std::vector<std::vector<double>> times2(Bs, std::vector<double>(M));
    std::vector<std::vector<double>> times3(Bs, std::vector<double>(M));


    for (size_t m = 0; m < M; m++) {

        printf("Working with matrices of size %d\n", matrixSize[m]);
        printf("---------------------------------------------\n");

        size_t N = matrixSize[m];
        std::vector<double> A(N * N), B(N * N), C(N * N);

        // TODO: Question 2c: Initialize matrices
        //       store A and B as row major and C as column major
        // Store A as row major
        for (size_t i=0; i<N; ++i)
            for (size_t j=0; j<N; ++j)
                A[i*N+j] = 2*i+j;

        // Store B as row major
        for (size_t i=0; i<N; ++i)
            for (size_t j=0; j<N; ++j)
                B[i*N+j] = 2*i+j;

        // Store C as column major
        for (size_t i=0; i<N; ++i)
            for (size_t j=0; j<N; ++j)
                C[j*N+i] = 2*i+j;

        printf("Start C=A*B (non optimized).\n");
        times1[m] = benchmark_AB(A, B, 1, 0, Ns);

        printf("---------------------------------------------\n");

        for (size_t b = 0; b < Bs; b++) {
            printf("Start C=A*B (optimized, row major, block size=%zu).\n", blockSize[b]);
            times2[b][m] = benchmark_AB(A, B, 2, blockSize[b], Ns);
        }

        printf("---------------------------------------------\n");

        for (size_t b = 0; b < Bs; b++) {
            printf("Start C=A*B (optimized, column major, block size=%zu).\n", blockSize[b]);
            times3[b][m] = benchmark_AB(A, C, 3, blockSize[b], Ns);
        }

        printf("==================================================\n");
    }


    FILE *fp = nullptr;
    fp = fopen("matrix_matrix_times.txt", "w");
    // write header to the file
    std::string header = "# N   unopt ";
    for (size_t b = 0; b < Bs; b++)
        header = header + "  br_" + std::to_string(blockSize[b]);
    for (size_t b = 0; b < Bs; b++)
        header = header + "  bc" + std::to_string(blockSize[b]);
    header = header + "\n";
    fprintf(fp, "%s", header.c_str());

    for (size_t m = 0; m < M; m++) {
        fprintf(fp, "%d %lf ", matrixSize[m], times1[m]);
        for (size_t b = 0; b < Bs; b++)
            fprintf(fp, "%lf ", times2[b][m]);
        for (size_t b = 0; b < Bs; b++)
            fprintf(fp, "%lf ", times3[b][m]);
        fprintf(fp, "\n");
    }
    fclose(fp);

  return 0;
}
