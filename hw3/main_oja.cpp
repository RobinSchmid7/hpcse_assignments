#include <random>
#include <algorithm>

#include <vector>
#include "code/utils.h"
#include "code/perceptron.h"

#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>

#include <omp.h>

int main (int argc, char** argv)
{

  ///////////////////////////////////////////////////////////////////////////
  // The implementation of Oja's and Strange's rule are sensitive to
  // initialization of the weights and the gradient might explode due to
  // overflows (nans in the code).
  ///////////////////////////////////////////////////////////////////////////


  // DATA PARAMETERS
  int D = 2;  // Data dimension
  int N = 1024; // Number of training samples
  int num_comp = 2; // Number of principal components
  std::string data_name = "2D"; // Data path
  std::string scaler = "center"; // Scaler type
  std::string weight_init = "normal"; // "normal" or "allsame"
  std::string method_name = "OJA"; // Method
  std::string data_path = "./data/"+data_name+"_dataset.txt"; // Data path
  // TRAINING PARAMETERS
  const int nepoch = 20000; // Number of epochs
  const double learn_rate = 1e-7; // Learning rate
  const double tolerance = 0.0; // 1e-18;
  const int batch_size = 1; // Batch-size
  const int check_every = 10; // Frequency of checking the convergence criterion

  // // DATA PARAMETERS
  // int D = 1850;  // Data dimension
  // int N = 1280; // Number of training samples
  // int num_comp = 5; // Number of principal components
  // std::string data_name = "faces"; // Data path
  // std::string scaler = "standard"; // Scaler type
  // std::string weight_init = "allsame"; // "normal" or "allsame"
  // std::string method_name = "OJA"; // Method
  // std::string data_path = "./data/"+data_name+"_dataset.txt"; // Data path
  // // TRAINING PARAMETERS
  // const int nepoch = 100;         // Number of epochs
  // const double learn_rate = 1e-6; // Learning rate
  // const double tolerance = 0.0; // Convergence criterion (tolerance)
  // const int batch_size = 32;      // Batch-size
  // const int check_every = 1; // Frequency of checking the convergence criterion


  double time_start;
  double time_end;

  /////////////////////////////////////////////////////////////////////////
  // PERCEPTRON OJA'S RULE IMPLEMENTATION
  // 1. Implement the forward method (forward pass) of the perceptron
  // 2. Implement Oja's rule
  // 3. Implement sanger's rule
  // 4. Compute the learned eigenvalues (read form std at the output)
  // 5. Reconstructed the data from the compressed representation
  /////////////////////////////////////////////////////////////////////////

  std::cout << "Loading dataset from: " << data_path << std::endl;
  double *data = utils::loadDataset(data_path, N, D);

  /////////////////////////////////////////////////////////////////////////
  // Transpose data to column major
  double *data_T = new double[N * D];
  utils::transposeData(data_T, data, N, D);

  /////////////////////////////////////////////////////////////////////////
  // Compute mean and standard deviation of your data
  time_start = omp_get_wtime();
  double *data_mean = new double[D];
  double *data_std = new double[D];
  utils::computeMean(data_mean, data_T, N, D);
  utils::computeStd(data_std, data_mean, data_T, N, D);
  utils::writeRowMajorMatrixToFile("./results/data/"+data_name+method_name+"_mean.txt", data_mean, 1, D);
  utils::writeRowMajorMatrixToFile("./results/data/"+data_name+method_name+"_std.txt", data_std, 1, D);

  // Time tracking
  time_end = omp_get_wtime();
  std::cout << "MEAN/STD TIME = " << time_end-time_start << " seconds\n";

  /////////////////////////////////////////////////////////////////////////
  // Normalize the data
  if(scaler.compare("standard") == 0){
    utils::standardizeRowMajor(data, data_mean, data_std, N, D);
  }
  else if (scaler.compare("center") == 0){
    utils::centerDataRowMajor(data, data_mean, N, D);
  }

  utils::writeRowMajorMatrixToFile("./results/data/"+data_name+method_name+"_data_example_test_PANTELIS.txt", data, N, D);


  ///////////////////////////////////////////////////////////////////////////
  // Random number generator to shuffle the dataset.
  ///////////////////////////////////////////////////////////////////////////
  struct {
    std::mt19937 gen;
    inline size_t operator()(size_t n) {
      std::uniform_int_distribution<size_t> dist(0, n ? n-1 : 0);
      return dist(gen);
    }
  } generator;

  /////////////////////////////////////////////////////////////////////////
  // Construct the perceptron

  // Number of inputs
  const int nInputs = D;
  // Number of outputs (hidden neurons)
  const int nOutputs = num_comp;
  assert(N % batch_size == 0);
  const int num_batches = (int)N/batch_size;
  std::cout << "Number of batches: " << num_batches << std::endl;

  // Create Network:
  Perceptron perceptron(nInputs, nOutputs, N, weight_init);


  // Checking the norm of the total weight matrix initially |W|
  // W=[w_1, w_2, ... ]
  double norm_ = utils::computeArrayMatrixNorm(perceptron.weights, nOutputs*nInputs);
  std::cout << "INITIAL NORM: Weight norm |w|= " << norm_  <<  std::endl;

  // Plot the norm of the individual components w_1, w_2, ...
  utils::plotComponentNorms(perceptron.weights, num_comp, D);


  /////////////////////////////////////////////////////////////////////////
  // Construct input batch

  // Construct input batch
  double *batch_input = new double[batch_size * nInputs];
  for (int iepoch = 1; iepoch < nepoch; iepoch++)
  {
    // Sample randomly indexes from the data
    std::vector<int> sample_ids(N);
    // Fill array: 0, 1, ..., N-1
    std::iota(sample_ids.begin(), sample_ids.end(), 0);
    std::random_shuffle(sample_ids.begin(), sample_ids.end(), generator);

    // Iterate over all batches
    for (int batch_num=0; batch_num < num_batches; ++batch_num){
      // Construct batch
      for (int k=0; k<batch_size; ++k)
      {
        // Accumulate data indexes (given by sample_ids vector)
        int n_on_batch = sample_ids[batch_num * batch_size + k];
        for (int d=0; d<nInputs;++d)
        {
          batch_input[k * nInputs + d] = data[n_on_batch * nInputs + d];
        }
      }


      ////////////////////////////////////////////////////////////////
      // 1. The forward method needed to perform the gradient step is
      // missing. Implement it.
      // TODO: Implement the forward method inside the Perceptron class
      // TODO:

      // Gradient computed inside the perceptron for Hebb's rule
      perceptron.hebbsRuleGradient(batch_input, batch_size);

      ////////////////////////////////////////////////////////////////
      // 2. Implement Oja's rule

      // TODO:
      // perceptron.ojasRuleGradient(batch_input, batch_size);
      // :TODO


      ////////////////////////////////////////////////////////////////
      // 3. Implement Sanger's rule

      // TODO:
      // perceptron.sangersRuleGradient(batch_input, batch_size);
      // :TODO

      // Utilities to print the norm of the gradient for debugging
      // perceptron.printGradientNorm();

      // Gradient normalization (maybe useful)
      // perceptron.normalizeGradient();

      // perceptron.printGradientNorm();

      // Parameters update
      perceptron.updateParams(learn_rate);
    }

    ////////////////////////////////////////////////////////////////
    // Checking convergence criterion every check_every epochs
    if(iepoch % check_every ==0)
    {
      // Computing the difference in the weights
      double norm_diff = utils::computeArrayMatrixDifferenceNorm(perceptron.weights_prev, perceptron.weights, nOutputs*nInputs);
      // Computing the norm of the total weight matrix |W|
      double norm_ = utils::computeArrayMatrixNorm(perceptron.weights, nOutputs*nInputs);

      // Saving the learned components
      utils::writeColMajorMatrixToFile("./results/data/"+data_name+method_name+"_components.txt", perceptron.weights, nOutputs, nInputs);

      // Print statistics
      std::cout << "EPOCH: " << std::to_string(iepoch) << ", Weight norm |w|= " << norm_ << ", Difference |Dw|: " << norm_diff <<  std::endl;
      // Plot component norms
      utils::plotComponentNorms(perceptron.weights, num_comp, D);

      // If the weight change smaller than the tolerance, terminate
      if(norm_<tolerance){
        break;
      }
    }
  }

  // SAVING THE WEIGHTS
  utils::writeColMajorMatrixToFile("./results/data/"+data_name+method_name+"_components.txt", perceptron.weights, nOutputs, nInputs);

  ////////////////////////////////////////////////////////////////
  // 4. Compute the learned eigenvalues

  // TODO: Compute the learned eigenvalues
  perceptron.computeEigenvalues(data, N);
  // :TODO

  // Write eigenvalues to file
  utils::writeRowMajorMatrixToFile("./results/data/" + data_name + method_name + "_eig.txt", perceptron.eigenvalues, 1, nOutputs);

  // Compressed data at the output
  double *output = perceptron.forward(data, N);

  ////////////////////////////////////////////////////////////////
  // 5. Reconstruct the original data
  // TODO:
  double *data_rec = new double[N * D];
  utils::reconstructDatasetColMajor(data_rec, perceptron.weights, output, N, nInputs, nOutputs);
  // :TODO

  /////////////////////////////////////////////////////////////////////////
  // Un-Normalize the data
  if(scaler.compare("standard") == 0){
    utils::inverseStandarizeDatasetRowMajor(data_rec, data_mean, data_std, N, D);
  }
  else if (scaler.compare("center") == 0){
    utils::inverseCenterDatasetRowMajor(data_rec, data_mean, N, D);
  }

  utils::writeRowMajorMatrixToFile("./results/data/"+data_name+method_name+"_data_reconstructed.txt", data_rec, N, D);

  delete[] batch_input;
  delete[] data;
  delete[] data_T;
  delete[] data_mean;
  delete[] data_std;
  delete[] data_rec;

  return 0;
}
