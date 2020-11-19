/*Author: Robin Schmid, 16-927-725, schmirob@ethz.ch*/
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <chrono>
#include <math.h>
#include <string>


void transpose( std::vector<double> &A ){

  size_t N = sqrt(A.size());
  std::vector<double> AT(N*N);
  // TODO: Question 2b: Straightforward matrix transposition
  // Swap i, j for transposition, naive way with row major
    for (size_t i=0; i<N; ++i)
      for (size_t j=0; j<N; ++j)
          AT[i*N+j] = A[j*N+i];
}


void transpose_block( std::vector<double> &A, size_t blockSize ){
  size_t N = sqrt(A.size());
  std::vector<double> AT(N*N);
  size_t b = blockSize;
  // TODO: Question 2b: Block matrix transposition
  // Iterate through all blocks
  for (size_t i=0; i<N/b; ++i)
    for (size_t j=0; j<N/b; ++j)
      // Iterate through all elements within a block
      for (size_t bi=0; bi<b; ++bi)
        for (size_t bj=0; bj<b; ++bj)
          // Swap index of blocks and elements within a block
          AT[j*N*b+i*b+bi*N+bj] = A[i*N*b+j*b+bj*N+bi];
}




double benchmark_transpose( std::vector<double> &A, size_t mode, size_t blockSize, size_t Ns ){

  size_t N = sqrt(A.size());
  double times = 0;

  // TODO: Check that the matrix size is divided by the blockSize when mode==2
  if( mode==2 &&  N%blockSize!=0 ){
    printf("Error: the size of the matrix (%zu) should be divided by the blockSize variable (%zu).\n",N,blockSize);
    exit(1);
  }

  for( size_t i=0; i<Ns; i++){
    auto t1 = std::chrono::system_clock::now();
    // TODO: Question 2b: Call the function to be benchmarked
    if( mode==1 ){
        transpose(A);
    }
    else if( mode==2 ){
        transpose_block(A, blockSize);
    }
    auto t2 = std::chrono::system_clock::now();
    times += std::chrono::duration<double>(t2-t1).count();
  }
  printf("Done in total %9.4fs  --  average %9.4fs\n", times, times/Ns);

  return times/Ns;

}


int main( ){

    std::vector<int> matrixSize{1024, 2048, 4096};
    size_t M = matrixSize.size();

    std::vector<size_t> blockSize{2, 4, 8, 16, 32, 64, 128};
    size_t B = blockSize.size();

    size_t Ns = 2;

    std::vector<double> times1(M);
    std::vector<std::vector<double>> times2(B, std::vector<double>(M));


    // loop over matrix sizes
    for (size_t m = 0; m < M; m++) {

        printf("Working with a matrix of size %d\n", matrixSize[m]);

        size_t N = matrixSize[m];

        std::vector<double> A(N*N);

        // TODO: Initialize the matrix
        // Initialize similar to 2a) row major
        for (size_t i=0; i<N; ++i)
            for (size_t j=0; j<N; ++j)
                A[i*N+j] = i+j;

        printf("Start transposing (non optimized).\n");
        times1[m] = benchmark_transpose(A, 1, 0, Ns);

        // loop over block sizes
        for (size_t b = 0; b < B; b++) {

            printf("Start transposing (optimized, block size=%zu).\n", blockSize[b]);
            times2[b][m] = benchmark_transpose(A, 2, blockSize[b], Ns);
        }

        printf("==================================================\n");
    }


    // write results to a file
    FILE *fp = nullptr;
    fp = fopen("transpose_times.txt", "w");
    // write header to the file
    std::string header = "# N   time_unoptimized ";
    for (size_t b = 0; b < B; b++)
        header = header + "  block_" + std::to_string(blockSize[b]);
    header = header + "\n";
    fprintf(fp, "%s", header.c_str());
    for (size_t m = 0; m < M; m++) {
        fprintf(fp, "%d %lf ", matrixSize[m], times1[m]);
        for (size_t b = 0; b < B; b++)
            fprintf(fp, "%lf ", times2[b][m]);
        fprintf(fp, "\n");
    }
    fclose(fp);

  return 0;
}
