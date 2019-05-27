#include "cuda_runtime.h"
#include "device_launch_parameters.h"

//Riduzione di n in modulo p.
__device__ int mod_long_GPU(long long n, long long p) {
	long long v = n, x = 0;

	if (v >= p) {
		v = n % p;
	}
	else {
		if (v < 0) {
			x = n / p;
			v = n - (x*p);
			v += p;
		}
	}
	int r = v;
	return r;
}

//Scambio di due righe della matrice m.
__device__ void swap_rows_GPU(int *m, int row, int col, int j, int i) {

	int k;
	long long tmp;
	if (j != i) {
		for (k = 0;k<col;k++) {
			tmp = m[i*col + k];				//m[i][k];
			m[i*col + k] = m[j*col + k];	//m[i][k] = m[j][k];
			m[j*col + k] = tmp;				//m[j][k] = tmp;
		}
	}
}



//Riduzione di n in modulo p.
__device__ int mod_GPU(int n, int p) {
	int v = n, x = 0;

	if (v >= p) {
		v = n % p;
	}
	else {
		if (v < 0) {
			x = n / p;
			v = n - (x*p);
			v += p;
		}
	}
	return v;
}


//inverso moltiplicativo di n in modulo p (con p primo).
__device__ int invers_GPU(int n, int p) {
	int b0 = p, t, q;
	int x0 = 0, x1 = 1;
	if (p == 1) return 1;
	while (n > 1) {
		q = n / p;
		t = p, p = (n % p), n = t;
		t = x0, x0 = x1 - q * x0, x1 = t;
	}
	if (x1 < 0) x1 += b0;
	return x1;
}


// a + b mod p
//sommatoria di a e b in modulo p
__device__ int add_mod_GPU(int a, int b, int p) {
	return mod_GPU((a + b), p);
}

// a - b mod p
//sottrazione di a e b in modulo p
__device__ int sub_mod_GPU(int a, int b, int p) {
	long long aa, bb;
	aa = a;
	bb = b;
	return mod_long_GPU((aa - bb), p);
}

// a * b mod p
//prodotto di a e b in modulo p
__device__ int mul_mod_GPU(int a, int b, int p) {
	long long aa, bb;
	aa = a;
	bb = b;
	return mod_long_GPU((aa*bb), p);
}

__global__ void reset_pivot_col(int *matrix, int row, int col, int pivot_row, int pivot_col, int thread_height, int block_dim) {

	int start_row = (pivot_row + 1) + ((blockIdx.x * (thread_height*block_dim)) + (threadIdx.x * thread_height));
	int reached_row = (pivot_row + 1) + ((blockIdx.x * (thread_height*block_dim)) + ((threadIdx.x + 1) * thread_height));
	int iteration = thread_height;
	if (reached_row > row) {
		iteration = thread_height - (reached_row - row);
		if (iteration > thread_height) {
			iteration = 0;
		}
	}

	for (int i = 0; i<iteration; i++) {
		matrix[(start_row + i)*col + pivot_col] = 0;
	}
}

__global__ void swap_rows(int *matrix, int row, int col, int j, int i) {

	int tid = blockIdx.x * blockDim.x + threadIdx.x;
	if (tid >= col) {
		return;
	}
	int ii = i * col + tid;
	int jj = j * col + tid;
	int tmp = matrix[ii];
	matrix[ii] = matrix[jj];
	matrix[jj] = tmp;
}


