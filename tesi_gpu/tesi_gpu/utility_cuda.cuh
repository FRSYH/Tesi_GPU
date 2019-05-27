#include "cuda_runtime.h"
#include "device_launch_parameters.h"

__device__ int mod_long_GPU(long long n, long long p);

__device__ void swap_rows_GPU(int *m, int row, int col, int j, int i);

__device__ int mod_GPU(int n, int p);

__device__ int invers_GPU(int n, int p);

__device__ int add_mod_GPU(int a, int b, int p);

__device__ int sub_mod_GPU(int a, int b, int p);

__device__ int mul_mod_GPU(int a, int b, int p);

__global__ void reset_pivot_col(int *matrix, int row, int col, int pivot_row, int pivot_col, int thread_height, int block_dim);

__global__ void swap_rows(int *matrix, int row, int col, int j, int i);

