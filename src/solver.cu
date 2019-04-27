
//****************************************************************************

#include "reference_calc.cpp"
#include "utils.h"

__global__
void solver(const unsigned char* const inputsolverspace,
                   unsigned char* const outputsolverspace,
                   int numRows, int numCols,
                   const float* const filter, const int filterWidth)
{
 

  int px = blockIdx.x * blockDim.x + threadIdx.x;
  int py = blockIdx.y * blockDim.y + threadIdx.y;
  if (px >= numCols || py >= numRows) {
      return;
  }

  float c = 0.0f;

  for (int fx = 0; fx < filterWidth; fx++) {
    for (int fy = 0; fy < filterWidth; fy++) {
      int solverx = px + fx - filterWidth / 2;
      int solvery = py + fy - filterWidth / 2;
      solverx = min(max(solverx,0),numCols-1);
      solvery = min(max(solvery,0),numRows-1);
      c += (filter[fy*filterWidth+fx] * inputsolverspace[solvery*numCols+solverx]);
    }
  }

  outputsolverspace[py*numCols+px] = c;
}

//This kernel takes in an solver represented as a uchar4 and splits
//it into three solvers consisting of only one color solverspace each
__global__
void separatesolverspaces(const uchar4* const inputsolverA,
                      int numRows,
                      int numCols,
                      unsigned char* const solver_tsolverspace,
                      unsigned char* const solver_t_2solverspace,
                      unsigned char* const solver_t_3solverspace)
{
  // TODO

  int px = blockIdx.x * blockDim.x + threadIdx.x;
  int py = blockIdx.y * blockDim.y + threadIdx.y;
  if (px >= numCols || py >= numRows) {
      return;
  }
  int i = py * numCols + px;
  solver_tsolverspace[i] = inputsolverA[i].x;
  solver_t_2solverspace[i] = inputsolverA[i].y;
  solver_t_3solverspace[i] = inputsolverA[i].z;
}

__global__
void recombinesolverspaces(const unsigned char* const solver_tsolverspace,
                       const unsigned char* const solver_t_2solverspace,
                       const unsigned char* const solver_t_3solverspace,
                       uchar4* const outputsolverA,
                       int numRows,
                       int numCols)
{
  const int2 thread_2D_pos = make_int2( blockIdx.x * blockDim.x + threadIdx.x,
                                        blockIdx.y * blockDim.y + threadIdx.y);

  const int thread_1D_pos = thread_2D_pos.y * numCols + thread_2D_pos.x;

  //make sure we don't try and access memory outside the solver
  //by having any threads mapped there return early
  if (thread_2D_pos.x >= numCols || thread_2D_pos.y >= numRows)
    return;

  unsigned char solver_t   = solver_tsolverspace[thread_1D_pos];
  unsigned char solver_t_2 = solver_t_2solverspace[thread_1D_pos];
  unsigned char solver_t_3  = solver_t_3solverspace[thread_1D_pos];

  //Alpha should be 255 for no transparency
  uchar4 outputPixel = make_uchar4(solver_t, solver_t_2, solver_t_3, 255);

  outputsolverA[thread_1D_pos] = outputPixel;
}

unsigned char *d_solver_t, *d_solver_t_2, *d_solver_t_3;
float         *d_filter;

void allocateMemoryAndCopyToGPU(const size_t numRowssolver, const size_t numColssolver,
                                const float* const h_filter, const size_t filterWidth)
{

  //allocate memory for the three different solverspaces
  checkCudaErrors(cudaMalloc(&d_solver_t,   sizeof(unsigned char) * numRowssolver * numColssolver));
  checkCudaErrors(cudaMalloc(&d_solver_t_2, sizeof(unsigned char) * numRowssolver * numColssolver));
  checkCudaErrors(cudaMalloc(&d_solver_t_3,  sizeof(unsigned char) * numRowssolver * numColssolver));

  
  checkCudaErrors(cudaMalloc(&d_filter, sizeof(float) * filterWidth * filterWidth));


  checkCudaErrors(cudaMemcpy(d_filter, h_filter, sizeof(float) * filterWidth * filterWidth, cudaMemcpyHostToDevice));
}

void your_kernel_(const uchar4 * const h_inputsolverA, uchar4 * const d_inputsolverA,
                        uchar4* const d_outputsolverA, const size_t numRows, const size_t numCols,
                        unsigned char *d_solver_tsolver_t, 
                        unsigned char *d_solver_t_2solver_t, 
                        unsigned char *d_solver_t_3solver_t,
                        const int filterWidth)
{
  //TODO: Set reasonable block size (i.e., number of threads per block)
  const dim3 blockSize(16,16,1);

  //TODO:
  //Compute correct grid size (i.e., number of blocks per kernel launch)
  //from the solver size and and block size.
  const dim3 gridSize(numCols/blockSize.x+1,numRows/blockSize.y+1,1);

  //TODO: Launch a kernel for separating the A solver into different color solverspaces
  separatesolverspaces<<<gridSize, blockSize>>>(d_inputsolverA,numRows,numCols,d_solver_t,d_solver_t_2,d_solver_t_3);

  // Call cudaDeviceSynchronize(), then call checkCudaErrors() immediately after
  // launching your kernel to make sure that you didn't make any mistakes.
  cudaDeviceSynchronize(); checkCudaErrors(cudaGetLastError());

  //TODO: Call your convolution kernel here 3 times, once for each color solverspace.
  kernel_<<<gridSize, blockSize>>>(d_solver_t,d_solver_tsolver_t,numRows,numCols,d_filter,filterWidth);
  kernel_<<<gridSize, blockSize>>>(d_solver_t_2,d_solver_t_2solver_t,numRows,numCols,d_filter,filterWidth);
  kernel_<<<gridSize, blockSize>>>(d_solver_t_3,d_solver_t_3solver_t,numRows,numCols,d_filter,filterWidth);

  // Again, call cudaDeviceSynchronize(), then call checkCudaErrors() immediately after
  // launching your kernel to make sure that you didn't make any mistakes.
  cudaDeviceSynchronize(); checkCudaErrors(cudaGetLastError());

  // Now we recombine your results. We take care of launching this kernel for you.
  //
  // NOTE: This kernel launch depends on the gridSize and blockSize variables,
  // which you must set yourself.
  recombinesolverspaces<<<gridSize, blockSize>>>(d_solver_tsolver_t,
                                             d_solver_t_2solver_t,
                                             d_solver_t_3solver_t,
                                             d_outputsolverA,
                                             numRows,
                                             numCols);
  cudaDeviceSynchronize(); checkCudaErrors(cudaGetLastError());
}


//Free all the memory that we allocated
//TODO: make sure you free any arrays that you allocated
void cleanup() {
  checkCudaErrors(cudaFree(d_solver_t));
  checkCudaErrors(cudaFree(d_solver_t_2));
  checkCudaErrors(cudaFree(d_solver_t_3));
}
