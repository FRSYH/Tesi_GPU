#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <time.h>
#include <stdbool.h>
#include <time.h>
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include "linalg.h"
#include "matrix.h"
#include "scan.h"
#include "utility.h"
#include "utility_cuda.cuh"

// compilazione nvcc gm.cu -o gm -w -Xcompiler " -openmp"
// nvcc gm.cu -o gm -w -Xcompiler " -openmp" -gencode arch=compute_61,code=sm_61 -lcudadevrt -rdc=true

__device__ int next_pivot_row = 0;

//dichiarazione variabili globali
int max_degree = 0;
int module = 0;

//----------------------------------------------------------------------------------------------------------------
#define gpuErrchk(ans) { gpuAssert((ans), __FILE__, __LINE__); }
inline void gpuAssert(cudaError_t code, const char *file, int line, bool abort = true)
{
	if (code != cudaSuccess)
	{
		fprintf(stderr, "GPUassert: %s %s %d\n", cudaGetErrorString(code), file, line);
		if (abort) exit(code);
	}
}

__global__ void find_pivot(int *matrix, int row, int col, int r, int pivot_col) {

	int tid = blockIdx.x * blockDim.x + threadIdx.x;
	int thread_row = r + tid;
	if (thread_row >= row)
		return;
	if (matrix[thread_row*col + pivot_col] != 0) {
		atomicMin(&next_pivot_row, thread_row);
	}
}

__global__ void submatrix_reduction_by_cell(int *matrix, int row, int col, int module, int inv, int pivot_col, int pivot_row) {
	
	int starting_row = pivot_row + 1;
	int idx = blockIdx.x * blockDim.x + threadIdx.x;
	int idy = blockIdx.y * blockDim.y + threadIdx.y + starting_row;

	if (idx < pivot_col  && idy < row && idy > pivot_row) {		//fermo i thread prima di pivot_colonna per impedire di sovrascrivere il dato necessario per s
		int div = matrix[idy*col + pivot_col];
		if (div != 0) {
			int s = mul_mod_GPU(inv, div, module);
			int a = mul_mod_GPU(s, matrix[pivot_row*col + idx], module);
			matrix[idy*col + idx] = sub_mod_GPU(matrix[idy*col + idx], a, module);
		}
	}
}


__global__ void gaussian_reduction_by_cell(int *matrix, int row, int col, int module) {

	int pivot_row = 0, r = 0, rows_found = 0;
	int inv;

	for (int pivot_col = col - 1; pivot_col >= 0; pivot_col--) {
		r = rows_found;
		while (r < row && matrix[r*col + pivot_col] == 0) {   //m[r][pivot_colonna]
			r++;
		}
		// ho trovato la prima riga con elemento non nullo in posizione r e pivot_colonna oppure non esiste nessuna riga con elemento non nullo in posizione pivot_colonna

		if (r < row) { //significa che ho trovato un valore non nullo
			if (r != rows_found) {
				swap_rows_GPU(matrix, row, col, rows_found, r); //sposto la riga appena trovata nella posizone corretta
			}
			pivot_row = rows_found;
			rows_found++;

			inv = invers_GPU(matrix[pivot_row*col + pivot_col], module);		//inverso dell´ elemento in m[r][pivot_colonna]	

																					//kernel per riduzione celle
			int block_dim = 16;
			dim3 threads(block_dim, block_dim, 1);
			int number_of_rows = row - rows_found;
			int grid_y = number_of_rows / block_dim + 1;
			int grid_x = col / block_dim + 1;
			dim3 blocks(grid_x, grid_y, 1);

			submatrix_reduction_by_cell<<<blocks, threads>>>(matrix, row, col, module, inv, pivot_col, pivot_row);
			cudaDeviceSynchronize();

			//necessario azzerare tutta la colonna (pivot_colonna)
			for (int x = pivot_row + 1; x < row; x++) {
				matrix[x*col + pivot_col] = 0;
			}
		}
	}
}

__global__ void submatrix_reduction_by_row(int *matrix, int row, int col, int module, int start, int pivot_col, int inv, int pivot_row, int cell_per_thread) {

	extern __shared__ int smem[];
	if ((threadIdx.x * cell_per_thread) <= pivot_col) {
		int row_offset = pivot_row * col;
		int thread_offset = threadIdx.x * cell_per_thread;
		//allocazione della smem con la riga di pivot, ogni thread copia una porzione di riga pari a "cell_per_thread".
		for (int i = 0; i<cell_per_thread; i++) {
			if (thread_offset + i <= pivot_col) {
				smem[thread_offset + i] = matrix[row_offset + thread_offset + i];
			}
		}
	}

	__syncthreads();

	int x = 0, y = 0;
	int row_index = (pivot_row + 1) + (blockDim.x * blockIdx.x + threadIdx.x);
	if (row_index >= start && row_index < row) {

		int row_linear_index = row_index * col + pivot_col;
		if (matrix[row_linear_index] != 0) {
			y = mul_mod_GPU(inv, matrix[row_linear_index], module);
			for (int k = 0; k < pivot_col + 1; k++) {
				//a = mul_mod_GPU(s,matrix[pivot_riga*col+k],module);
				x = mul_mod_GPU(y, smem[k], module);
				matrix[row_index*col + k] = sub_mod_GPU(matrix[row_index*col + k], x, module);
			}
		}
	}
}

__global__ void gaussian_reduction_by_row(int *matrix, int row, int col, int module, int dim) {

	int pivot_row = 0, r = 0, rows_found = 0;
	int inv;

	for (int pivot_col = col - 1; pivot_col >= 0; pivot_col--) {
		r = rows_found;
		while (r < row && matrix[r*col + pivot_col] == 0) {   //m[r][pivot_colonna]
			r++;

		}
		// ho trovato la prima riga con elemento non nullo in posizione r e pivot_colonna oppure non esiste nessuna riga con elemento non nullo in posizione pivot_colonna

		if (r < row) { //significa che ho trovato un valore non nullo
			if (r != rows_found) {
				swap_rows_GPU(matrix, row, col, rows_found, r); //sposto la riga appena trovata nella posizone corretta
			}
			pivot_row = rows_found;
			rows_found++;

			inv = invers_GPU(matrix[pivot_row*col + pivot_col], module);		//inverso dell´ elemento in m[r][pivot_colonna]	

			int block_dim = 1024;
			//kernel per riduzione righe	
			int numero_righe = row - rows_found;
			int t = (numero_righe < block_dim ? numero_righe : block_dim);
			int b = 1;
			if (t == block_dim && numero_righe != block_dim) {
				b = numero_righe / block_dim + 1;
			}

			dim3 threads(t);
			dim3 blocks(b);

			int pivot_length = pivot_col + 1;
			int cell_per_thread = (t >= pivot_length) ? 1 : (pivot_length / t) + 1;
			int shared_mem = pivot_length * sizeof(int);

			submatrix_reduction_by_row<<<blocks, threads, shared_mem >>>(matrix, row, col, module, rows_found, pivot_col, inv, pivot_row, cell_per_thread);
			cudaDeviceSynchronize();

		}
	}
}

__global__ void submatrix_reduction_by_block(int *matrix, int row, int col, int module, int pivot_col, int inv, int pivot_row, int thread_height, int block_dim) {

	extern __shared__ int smem[];

	int *smem_pivot_row = (int*)smem;
	int *smem_pivot_col = (int*)&smem_pivot_row[block_dim];

	int x = 0, y = 0, interation = 0;
	int col_index = blockIdx.x * blockDim.x + threadIdx.x;	//indice della colonna della matrice originale per il thread corrente
	//-------------
	//inizzializzazione smem per pivot riga
	smem_pivot_row[threadIdx.x] = matrix[pivot_row * col + col_index];	//ogni thread copia un solo elemento nella riga in shared, un thread per cella di riga
	//------------
	//inizializzazione smem per pivot colonna
	//calcolo del numero di celle (colonna_pivot) che ogni thred deve copiare
	
	int cell_to_copy = 1;
	if (thread_height > blockDim.x) {
		cell_to_copy = thread_height / blockDim.x + 1;
	}

	int base_row = (pivot_row + 1) + blockIdx.y * thread_height;
	int index = 0;
	//copia della porzione di colonna in smem
	for (int i = 0; i<cell_to_copy; i++) {
		index = (threadIdx.x * cell_to_copy) + i;
		if (base_row + index < row && index < thread_height) {
			smem_pivot_col[index] = matrix[(base_row + index) * col + pivot_col];
		}
	}
	//sincronizza tutti i thread del blocco in modo tale che la smem sia consistente
	__syncthreads();

	if (col_index < pivot_col) {
		//calcolo del numero di righe sulle quali deve iterare il thread, caso in cui la dimensione della matrice non collima con thread_height
		int reached_row = (pivot_row + 1) + ((blockIdx.y + 1) * thread_height); //riga raggiunta dal thread corrente
		if (reached_row > row) {
			interation = thread_height - (reached_row - row);	//dimensione non collima
		}
		else {
			interation = thread_height;	//caso normale
		}

		int row_offset = (pivot_row + 1) + (blockIdx.y * thread_height);

		for (int i = 0; i<interation; i++) {
			int pivot_element = smem_pivot_col[i];
			if (pivot_element != 0) {
				y = mul_mod_GPU(inv, pivot_element, module);		//tutti i thread sulla stessa riga calcolano lo stesso risultato
				x = mul_mod_GPU(y, smem_pivot_row[threadIdx.x], module);
				matrix[row_offset * col + col_index] = sub_mod_GPU(matrix[row_offset * col + col_index], x, module);
			}
			row_offset++;
		}
	}
}


__global__ void gaussian_reduction_by_block(int *matrix, int row, int col, int module, int dim) {

	int pivot_row = 0, r = 0, rows_found = 0;
	int inv;
	int block_dim = 0;
	int threads_per_block = 0;
	int block_x_axis, block_y_axis = 0;

	for (int pivot_col = col - 1; pivot_col >= 0; pivot_col--) {
		r = rows_found;
		///////////////////////////FIND PIVOT///////////////////////////////////////////////
		block_dim = 256;
		int row_to_check = row - rows_found;
		threads_per_block = (row_to_check < block_dim ? row_to_check : block_dim);
		dim3 t_find(threads_per_block);

		if (threads_per_block == block_dim && row_to_check != block_dim) {
			block_x_axis = (row_to_check / block_dim) + 1;
		}
		else {
			block_x_axis = 1;
		}
		dim3 b_find(block_x_axis);
		next_pivot_row = row;
		find_pivot << <b_find, t_find >> >(matrix, row, col, r, pivot_col);
		cudaDeviceSynchronize();
		r = next_pivot_row;
		/////////////////////////////////////////////////////////////////////////////////
		if (r < row) {
			if (r != rows_found) {
				////////////////////////SWAP ROWS////////////////////////////////////////////////////////
				block_dim = 256;
				threads_per_block = (col < block_dim ? col : block_dim);
				dim3 t_swap(threads_per_block);

				if (threads_per_block == block_dim && col != block_dim) {
					block_x_axis = (col / block_dim) + 1;
				}
				else {
					block_x_axis = 1;
				}
				dim3 b_swap(block_x_axis);
				//sposto la riga appena trovata nella posizone corretta
				swap_rows << <b_swap, t_swap >> >(matrix, row, col, rows_found, r);
				cudaDeviceSynchronize();

				////////////////////////////////////////////////////////////////////////////////////////
			}
			pivot_row = rows_found;
			rows_found++;

			inv = invers_GPU(matrix[pivot_row*col + pivot_col], module);
			////////////////////////////////////////REDUCTION BY BLOCK////////////////////////////////////
			block_dim = 128;
			int col_to_reduce = pivot_col;
			threads_per_block = (col_to_reduce < block_dim ? col_to_reduce : block_dim);
			dim3 threads(threads_per_block);

			if (threads_per_block == block_dim && col_to_reduce != block_dim) {
				block_x_axis = (col_to_reduce / block_dim) + 1;
			}
			else {
				block_x_axis = 1;
			}

			int thread_height = 256;
			int row_to_reduce = row - rows_found;
			block_y_axis = (row_to_reduce / thread_height) + 1;

			dim3 blocks(block_x_axis, block_y_axis);

			int shared = (block_dim * sizeof(int)) + (thread_height * sizeof(int));
			submatrix_reduction_by_block << <blocks, threads, shared >> >(matrix, row, col, module, pivot_col, inv, pivot_row, thread_height, block_dim);
			cudaDeviceSynchronize();
			//////////////////////////////////////////////////////////////////////////////////////////

			///////////////////////////////RESET PIVOT COL////////////////////////////////////////
			thread_height = 100;
			block_dim = 32;
			row_to_reduce = row - pivot_row;
			threads_per_block = (row_to_reduce < thread_height ? 1 : block_dim);
			block_x_axis = (threads_per_block == block_dim && row_to_reduce != block_dim) ? (row_to_reduce / (thread_height*block_dim) + 1) : 1;
			dim3 t(threads_per_block);
			dim3 b(block_x_axis);

			reset_pivot_col << <b, t >> >(matrix, row, col, pivot_row, pivot_col, thread_height, block_dim);
			cudaDeviceSynchronize();
			//////////////////////////////////////////////////////////////////////////////////////
		}
	}
}

double gauss_CUDA(int *m, int row, int col, int module) {

	int matrix_length = row * col;
	int matrix_length_bytes = matrix_length * sizeof(int);
	clock_t start, end;
	double elapsed = 0.0;
	int *m_d;

	gpuErrchk(cudaMalloc((void **)&m_d, matrix_length_bytes));
	gpuErrchk(cudaMemcpy(m_d, m, matrix_length_bytes, cudaMemcpyHostToDevice));
	start = clock();
	gaussian_reduction_by_block<<<1, 1>>>(m_d, row, col, module, row*col);
	gpuErrchk(cudaDeviceSynchronize());
	gpuErrchk(cudaMemcpy(m, m_d, matrix_length_bytes, cudaMemcpyDeviceToHost));

	end = clock();
	elapsed = ((double)(end - start)) / CLOCKS_PER_SEC;
	gpuErrchk(cudaFree(m_d));
	return elapsed;
}

void execute_standard(int **matrix, int * row, int col, struct map map, int *degree, int **monomi, int numero_variabili, int n_loops, int expansion, FILE *output_file) {

	clock_t start, end;
	double elapsed;

	//creo l'array che conterrà i gradi dei vari round
	int **m_deg_array = (int **)malloc(sizeof(int*));
	m_deg_array[0] = (int *)calloc(max_degree + 1, sizeof(int));
	int n_round = 0;
	int *m_deg = m_deg_array[0];
	int missing_degree = max_degree;
	fprintf(output_file, "Inizio computazione, metodo standard\n");
	matrix_degree(*matrix, *row, col, m_deg, monomi, numero_variabili);

	int stop = 0;

	while (stop != 1) {
		n_round++;

		fprintf(output_file, "\n -Eseguo moltiplicazione, ");
		fflush(stdout);

		start = clock();

		//find missing degree to multiply matrix
		for (int i = max_degree; i>0; i--) {
			if (m_deg[i] == 0) {
				missing_degree = i;
				break;
			}
		}

		moltiplica_matrice(matrix, row, col, map, degree, monomi, numero_variabili, missing_degree, max_degree);

		end = clock();
		elapsed = ((double)(end - start)) / CLOCKS_PER_SEC;
		fprintf(output_file, "numero righe: %d     (%f sec)", *row, elapsed);

		fprintf(output_file, "\n -Eseguo Gauss, ");
		fflush(stdout);
		//start = clock();

		//applico la riduzione di Gauss
		elapsed = gauss_CUDA(*matrix, *row, col, module);
		//elimino le righe nulle della matrice
		eliminate_null_rows(matrix, row, col);

		//aggiungo all'array i gradi dell'attuale round
		//n_round+1 perchè salvo anche i gradi prima di inziare i round
		m_deg_array = (int **)realloc(m_deg_array, sizeof(int*)*(n_round + 1));
		m_deg_array[n_round] = (int *)calloc(max_degree + 1, sizeof(int));
		m_deg = m_deg_array[n_round];

		//end = clock();
		//elapsed =  ((double)(end - start)) / CLOCKS_PER_SEC;
		fprintf(output_file, "numero righe: %d               (%f sec)\n", *row, elapsed);

		matrix_degree(*matrix, *row, col, m_deg, monomi, numero_variabili);
		print_matrix_degree(m_deg, output_file, max_degree);

		if (target_degree(m_deg, max_degree) == 0)
			stop = 1;

	}
	for (int i = 0; i < n_round + 1; i++)
		free(m_deg_array[i]);
	free(m_deg_array);
}


int main(int argc, char *argv[]) {

	FILE *input_file = NULL, *output_file = NULL;
	//inizializzo flag a false

	for (int parsed = 1; parsed < argc; parsed++) {
		if (parsed < argc && !strcmp(argv[parsed], "--input")) {
			parsed++;
			input_file = fopen(argv[parsed], "r");
			if (!input_file) {
				perror("Errore nell'apertura del file di input");
				return (-1);
			}
		}
		else if (parsed < argc && !strcmp(argv[parsed], "--output")) {
			parsed++;
			output_file = fopen(argv[parsed], "w");
			if (!output_file) {
				perror("Errore nell'apertura del file di output");
				return (-1);
			}
		}
	}

	if (!input_file)
		input_file = stdin;
	if (!output_file)
		output_file = stdout;

	int row, col, numero_variabili, tipo_ordinamento;
	int *matrix;

	char *variabili;
	row = col = numero_variabili = 0;
	int(*ord) (void*, const void *, const void *);
	struct map smap;

	clock_t start, end;
	double elapsed = 0.0;
	start = clock();

	//alloca la matrice principale, legge da input: il modulo,massimo grado e numero variabili
	allocation(&matrix, &row, &col, &numero_variabili, &variabili, &tipo_ordinamento, &module, &max_degree, input_file);

	if (order(&ord, tipo_ordinamento) != 0) {
		fprintf(stderr, "Ordinamento insesistente!!!\n\nTERMINAZIONE PROGRAMMA");
		return 0;
	}

	int * degree = (int *)calloc(max_degree + 1, sizeof(int));
	int numero_monomi = col;
	int **monomi;

	//crea il vettore con tutti i possibili monomi avendo num_var varaibili e max_degree come massimo grado
	monomi = monomial_computation(numero_variabili, max_degree, numero_monomi);

	//ordina il vettore dei monomi secondo un determinato ordinamento, ordinamento intercambiabile
	qsort_s(monomi, numero_monomi, sizeof(int*), ord, &numero_variabili);

	//inizializzazione matrice (lettura dati input)
	
	if (parse(numero_variabili, variabili, matrix, row, monomi, col, module, ord, input_file) == -1) {
		fprintf(stderr, "Errore di input !!!\n\nTERMINAZIONE PROGRAMMA"); //se l'input è in formato scorrettro abort del programma
		return 0;
	}

	end = clock();
	elapsed = ((double)(end - start)) / CLOCKS_PER_SEC;
	fprintf(output_file, "\nInizializzazione in %f sec\n", elapsed);

	start = clock();

	setup_struct_map(&smap, monomi, numero_monomi, numero_variabili, max_degree, ord);

	end = clock();
	elapsed = ((double)(end - start)) / CLOCKS_PER_SEC;
	fprintf(output_file, "\nMappa creata in %f sec,   %d x %d \n\n", elapsed, col, col);

	//RISOLUZIONE PROBLEMA
	start = clock();

	//inizializzazione vettore dei gradi dei polinomi
	init_degree_vector(degree, numero_variabili, max_degree);

	int n_loops = 30, expansion = 1;
	//eseguo moltiplicazione e riduzione di Gauss finche non trovo soluzione
	execute_standard(&matrix, &row, col, smap, degree, monomi, numero_variabili, n_loops, expansion, output_file);

	end = clock();
	elapsed = ((double)(end - start)) / CLOCKS_PER_SEC;
	fprintf(output_file, "\nTarget raggiunto, soluzione trovata in %f sec\n\n", elapsed);

	//print_matrix(matrix, row, col, output_file);
	print_incognite(matrix, row, col, numero_variabili, monomi, output_file);

	free(matrix);
	free(degree);
	cudaDeviceReset();

	return 0;
}

