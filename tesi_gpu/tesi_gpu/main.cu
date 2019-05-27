#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <time.h>
#include <stdbool.h>
#include <time.h>
#include <cuda_runtime.h>
#include <device_launch_parameters.h>
#include "linalg.h"
#include "matrix.h"
#include "scan.h"
#include "utility.h"

// compilazione nvcc gm.cu -o gm -w -Xcompiler " -openmp"
// nvcc gm.cu -o gm -w -Xcompiler " -openmp" -gencode arch=compute_61,code=sm_61 -lcudadevrt -rdc=true

__device__ int next_pivot_row = 0;

//dichiarazione variabili globali
int max_degree = 0;
int module = 0;

//----------------------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------------------


/*restituisce un array contenente tutti i len monomi con n variabili e grado <= m
len è il numero di possibili monomi con n variabili e grado <= m
i monomi sono array di interi di lunghezza n dove il valore di ogni posizione rappresenta
il grado della variabile in quella posizione. Esempio: n=3, x^2*y*z = [2,1,1]
len viene passato come argomento per evitare di ricalcolarlo internamente
*/
int **monomial_computation(int n, int m, int len) {

	int **vet, *monomial;

	//alloco la memoria per l'array
	matrix_alloc_int(&vet, len, n);

	//strutture di supporto necessarie per il calcolo
	monomial = (int *)malloc(n * sizeof(int));
	int pos = 0;

	//il calcolo è fatto dalla funzione ricorsiva correttemente parametrizzata
	monomial_computation_rec(n, m, vet, 0, monomial, &pos);

	free(monomial);

	return vet;
}



int init_matrix(int *matrix, int row, int col, int **vet_grd, char *variabili, int num_var, int(*ord) (void*, const void *, const void *), FILE *input_file) {
	//Inizializza la matrice principale (dei coefficienti) con i coefficienti dei polinomi forniti come input.
	return parse(num_var, variabili, matrix, row, vet_grd, col, module, ord, input_file);
}



void setup_struct_map(struct map *map, int **monomi, int len, int n, int m, int(*compar) (void*, const void *, const void *)) {

	int sum, index = len;

	//	inizializzo la struttura map, la mappa ha len righe.
	map->len = len;
	map->row = (map_row *)malloc(map->len * sizeof(struct map_row));

	//per ogni monomio in vet
	int row, col, i, v;
	for (row = 0; row < len; row++) {
		index = 0;
		//dichiarati dentro per la parallelizzazione
		int *temp = (int *)malloc(n * sizeof(int));
		int *save = (int *)calloc(len, sizeof(int));
		//provo a moltiplicarlo con ogni monomio in vet
		for (col = 0; col < len; col++) {
			sum = 0;
			//eseguo il prodotto (sum è la somma dei gradi)
			for (v = 0; v < n; v++) {
				temp[v] = monomi[row][v] + monomi[col][v];
				sum += temp[v];
			}
			//se il grado del prodotto > grado massimo tutti i restanti prodotti
			//su quella riga sono > grado massimo
			if (sum > m) {

				//	a questo punto col è l'indice del primo elemento della mappa che non è possibile rappresentare, quindi la riga row ha solo col numero di celle e non len come prima.
				index = col;
				break;
			}
			//altrimenti cerco il prodotto in vet e metto l'indice in save
			else {
				save[col] = int((int **)(bsearch_r((void *)&temp, (void *)monomi, len, (sizeof(int*)), compar, &n)) - monomi);
			}
		}

		//	terminato il ciclo sulle colonne posso inizializzare la struttura perchè conosco tutti gli elementi da inserire	
		//  la riga attuale ha esattamente index elementi diversi da -1, quindi la riga avrà lunghezza pari a index precedentemente calcolato
		//  alloco la riga con un array da index elementi

		map->row[row].len = index;
		map->row[row].col = (int *)malloc(map->row[row].len * sizeof(int));
		//	a questo map devo copiare gli elementi generati dento alla struttura

		for (i = 0; i<map->row[row].len; i++)
			map->row[row].col[i] = save[i];

		free(temp);
		free(save);
	}
}

void init_degree_vector(int *degree, int num_var) {
	//inizializza il vettore degree con il numero di monomi di grado i-esimo <= del grado massimo
	int i, c;
	for (i = 0; i<max_degree + 1; i++) {
		c = combination(num_var, i);
		degree[i] = c;
	}
}

int grado_monomio(int posizione, int **vet, int num_var) {
	//Calcola il grado del monomio a partire dalla posizione occupata nel vettore (ordinato) delle posizioni rispetto l'ordinamento scelto.
	//(la posizione occupata deve essere corretta).
	int i, grado;
	grado = 0;
	for (i = 0; i<num_var; i++) {
		grado += vet[posizione][i];
	}
	return grado;
}

void matrix_degree(int *m, int row, int col, int *m_deg, int **vet, int num_var) {
	//m_deg è un vettore che ha lunghezza pari al grado massimo.
	//la funzione calcola i gradi dei polinomi presenti nella matrice.
	//Ogni cella del vettore m_deg rappresenta un grado, se esso compare nella matrice allora viene impostato a 1 o altrimenti.

	int i, j, last, grado, linear_index = 0;
	for (i = 0; i<row; i++) {
		for (j = col - 1; j>0; j--) {
			linear_index = i * col + j;
			if (m[linear_index] != 0) {
				last = j;           //posizione dell'ultimo coefficiente della riga
				break;
			}
		}
		grado = grado_monomio(last, vet, num_var);
		m_deg[grado] = 1;
	}
}


void moltiplica_matrice(int **m, int *row, int col, struct map map, int * degree, int **vet, int num_var, int expansion_degree) {

	int riga;
	int grado_massimo_riga, grado_massimo_monomio, i, j, last, new_row = 0;
	last = -1;
	int linear_index = 0;
	long long total_dim = 0;
	int *last_index = (int*)calloc(*row, sizeof(int));
	int *numero_polinomi = (int*)calloc(*row, sizeof(int));
	int numero_nuovi_polinomi = 0;

	for (riga = 0; riga<*row; riga++) {
		for (i = col - 1; i>0; i--) {
			linear_index = riga * col + i;
			if ((*m)[linear_index] != 0) {  //(*m)[riga][i] != 0
				last = i;
				break;
			}
		}
		//risalgo al grado del monomio appena trovato
		//scorro la lista delle posizioni di inizio dei monomi con lo stesso grado

		last_index[riga] = last;

		if (last != -1) {

			grado_massimo_riga = grado_monomio(last, vet, num_var);

			//calcolo il grado massimo che deve avere il monomio per cui moltiplicare
			grado_massimo_monomio = max_degree - grado_massimo_riga;
			// a questo punto conosco per quanti monomi devo moltiplicare e quindi
			// conosco il numero di righe che devo aggiungere alla matrice
			if (expansion_degree != 0) {
				if (grado_massimo_monomio > expansion_degree) {
					grado_massimo_monomio = expansion_degree;
				}
			}

			for (i = 1; i<(grado_massimo_monomio + 1); i++) {
				new_row += degree[i];
				numero_nuovi_polinomi += degree[i];
			}
			numero_polinomi[riga] = numero_nuovi_polinomi;
			numero_nuovi_polinomi = 0;
		}
	}

	//--------------------------------------------------------------

	//printf("nuove righe %d, totale righe %d", new_row, (*row+new_row) );

	//ridimensionamento della matrice
	total_dim = (*row * col) + (new_row * col);
	*m = (int *)realloc(*m, total_dim * sizeof(int));
	//azzeramento delle nuove righe
	for (i = *row; i<*row + new_row; i++) {
		for (j = 0; j<col; j++) {
			(*m)[i*col + j] = 0;
		}
	}

	int len = *row;
	for (riga = 0; riga<len; riga++) {
		if (last_index[riga] != -1) {
			for (i = 1; i<(numero_polinomi[riga] + 1); i++) {     								//scorre tutti i monomi per i quali posso moltiplicare
				for (j = 0; j<(last_index[riga] + 1); j++) {     								//scorre fino all'ultimo elemento della riga
																								//(*m)[*row][ map.row[i].col[j] ] = (*m)[riga][j];  				//shift nella posizione corretta indicata dalla mappa
					linear_index = *row * col + map.row[i].col[j];
					(*m)[linear_index] = (*m)[riga*col + j];
				}
				*row = *row + 1;											//aumento del conteggio delle righe
			}
		}
	}

	free(last_index);
	free(numero_polinomi);
}

void moltiplica_riga_forn(int **m, int *row, int col, int riga, struct map map, int * degree, int **vet, int num_var, int stop_degree) {

	int grado_massimo_riga, grado_massimo_monomio, i, j, last, new_row;
	last = -1;
	int linear_index = 0;
	long long total_dim = 0;
	//cerco la posizione dell'ultimo coefficiente non nullo del polinomio rappresentato nella riga.
	for (i = col - 1; i>0; i--) {
		linear_index = riga * col + i;
		if ((*m)[linear_index] != 0) {  //(*m)[riga][i] != 0
			last = i;
			break;
		}
	}
	//risalgo al grado del monomio appena trovato
	//scorro la lista delle posizioni di inizio dei monomi con lo stesso grado
	if (last != -1) {

		grado_massimo_riga = grado_monomio(last, vet, num_var);

		//calcolo il grado massimo che deve avere il monomio per cui moltiplicare
		grado_massimo_monomio = max_degree - grado_massimo_riga;
		// a questo punto conosco per quanti monomi devo moltiplicare e quindi
		// conosco il numero di righe che devo aggiungere alla matrice
		new_row = 0;

		for (i = 1; i<(grado_massimo_monomio + 1); i++) {
			new_row += degree[i];
		}

		total_dim = (*row * col) + (new_row * col);
		*m = (int *)realloc(*m, total_dim * sizeof(int));
		//azzeramento delle nuove righe
		for (i = *row; i<*row + new_row; i++) {
			for (j = 0; j<col; j++) {
				(*m)[i*col + j] = 0;
			}
		}

		for (i = 1; i<(new_row + 1); i++) {     								//scorre tutti i monomi per i quali posso moltiplicare
			for (j = 0; j<(last + 1); j++) {     								//scorre fino all'ultimo elemento della riga
																				//(*m)[*row][ map.row[i].col[j] ] = (*m)[riga][j];  				//shift nella posizione corretta indicata dalla mappa
				linear_index = *row * col + map.row[i].col[j];
				(*m)[linear_index] = (*m)[riga*col + j];
			}
			*row = *row + 1;											//aumento del conteggio delle righe
		}
	}

}

int target_degree(int *v) {
	//Controlla se il vettore v rappresenta la condizione di terminazione con gradi completi {1,2,3,...,max_degree}
	//Se la condizione è soddisfatta return 0 altrimenti -1

	int i, flag;
	flag = 0;
	for (i = 1; i<max_degree + 1; i++) {
		if (v[i] != 1) {
			flag = -1;
			break;
		}
	}
	return flag;
}

void print_incognite(int *m, int row, int col, int num_var, int **vet, FILE *output_file) {

	int grado, last;

	for (int r = row - (num_var + 1); r<row; r++) {

		//cerca la posizione dell'ulitmo elemento non nullo della riga r
		for (int i = col - 1; i >= 0; i--) {
			if (m[r*col + i] != 0) { //m[r][i] != 0
				last = i;
				break;
			}
		}
		//calcola il grado della riga r
		grado = grado_monomio(last, vet, num_var);
		//se il grado della riga r è 1 allora stampa tutta la riga della matrice
		if (grado == 1) {
			for (int j = 0; j<last + 1; j++) {
				fprintf(output_file, "%d ", m[r*col + j]); //m[r][j]
			}
			fprintf(output_file, "\n\n");
		}
	}
	fprintf(output_file, "\n");
}


//Scambio di due righe della matrice m.
__device__ void swap_rows_GPU(int *m, int row, int col, int j, int i) {

	int k;
	long long tmp;
	if (j != i) {
		for (k = 0;k<col;k++) {
			tmp = m[i*col + k];			//m[i][k];
			m[i*col + k] = m[j*col + k];	//m[i][k] = m[j][k];
			m[j*col + k] = tmp;			//m[j][k] = tmp;
		}
	}
}

//n mod p 
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


//n mod p 
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



#define gpuErrchk(ans) { gpuAssert((ans), __FILE__, __LINE__); }
inline void gpuAssert(cudaError_t code, const char *file, int line, bool abort = true)
{
	if (code != cudaSuccess)
	{
		fprintf(stderr, "GPUassert: %s %s %d\n", cudaGetErrorString(code), file, line);
		if (abort) exit(code);
	}
}



__global__ void kernel_riduzione_blocco(int *matrix, int row, int col, int module, int pivot_col, int inv, int pivot_row, int thread_height, int block_dim) {

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



__global__ void kernel_riduzione_riga(int *matrix, int row, int col, int module, int start, int pivot_col, int inv, int pivot_row, int cell_per_thread) {


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


__global__ void kernel_riduzione_cella(int *matrix, int row, int col, int module, int inv, int pivot_col, int pivot_row) {
	
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


__global__ void gauss_kernel_celle(int *matrix, int row, int col, int module) {

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

			kernel_riduzione_cella<<<blocks, threads>>>(matrix, row, col, module, inv, pivot_col, pivot_row);
			cudaDeviceSynchronize();

			//necessario azzerare tutta la colonna (pivot_colonna)
			for (int x = pivot_row + 1; x < row; x++) {
				matrix[x*col + pivot_col] = 0;
			}

		}
	}
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

__global__ void find_pivot(int *matrix, int row, int col, int r, int pivot_col) {
	/*
	while( r < row && matrix[r*col+pivot_colonna] == 0 ){
	r++;
	}
	pointer_r = r;
	*/

	int tid = blockIdx.x * blockDim.x + threadIdx.x;
	int thread_row = r + tid;
	if (thread_row >= row)
		return;
	if (matrix[thread_row*col + pivot_col] != 0) {
		atomicMin(&next_pivot_row, thread_row);
	}
}



__global__ void gauss_kernel_blocco(int *matrix, int row, int col, int module, int dim) {

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
		find_pivot<<<b_find, t_find>>>(matrix, row, col, r, pivot_col);
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
				swap_rows<<<b_swap, t_swap>>>(matrix, row, col, rows_found, r);
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
			kernel_riduzione_blocco<<<blocks, threads, shared>>>(matrix, row, col, module, pivot_col, inv, pivot_row, thread_height, block_dim);
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

			reset_pivot_col<<<b, t >>>(matrix, row, col, pivot_row, pivot_col, thread_height, block_dim);
			cudaDeviceSynchronize();
			//////////////////////////////////////////////////////////////////////////////////////
		}
	}
}

__global__ void gauss_kernel_righe(int *matrix, int row, int col, int module, int dim) {

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

			kernel_riduzione_riga<<<blocks, threads, shared_mem >>>(matrix, row, col, module, rows_found, pivot_col, inv, pivot_row, cell_per_thread);
			cudaDeviceSynchronize();

		}
	}
}

double gauss_GPU(int *m, int row, int col, int module) {

	int matrix_length = row * col;
	int matrix_length_bytes = matrix_length * sizeof(int);
	clock_t start, end;
	double elapsed = 0.0;
	int *m_d;

	gpuErrchk(cudaMalloc((void **)&m_d, matrix_length_bytes));
	gpuErrchk(cudaMemcpy(m_d, m, matrix_length_bytes, cudaMemcpyHostToDevice));

	start = clock();

	gauss_kernel_blocco<<<1, 1>>>(m_d, row, col, module, row*col);
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

		moltiplica_matrice(matrix, row, col, map, degree, monomi, numero_variabili, missing_degree);

		end = clock();
		elapsed = ((double)(end - start)) / CLOCKS_PER_SEC;
		fprintf(output_file, "numero righe: %d     (%f sec)", *row, elapsed);

		fprintf(output_file, "\n -Eseguo Gauss, ");
		fflush(stdout);
		//start = clock();

		//applico la riduzione di Gauss
		elapsed = gauss_GPU(*matrix, *row, col, module);
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

		if (target_degree(m_deg) == 0)
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

	/*
	int peak_clk = 1;
	cudaError_t err = cudaDeviceGetAttribute(&peak_clk, cudaDevAttrClockRate, 0);
	if (err != cudaSuccess) {printf("cuda err: %d at line %d\n", (int)err, __LINE__); return 1;}
	printf("peak clock rate: %dkHz", peak_clk);
	*/

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
	if (init_matrix(matrix, row, col, monomi, variabili, numero_variabili, ord, input_file) == -1) {
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
	init_degree_vector(degree, numero_variabili);

	int n_loops = 30, expansion = 1;
	//eseguo moltiplicazione e riduzione di Gauss finche non trovo soluzione
	execute_standard(&matrix, &row, col, smap, degree, monomi, numero_variabili, n_loops, expansion, output_file);

	end = clock();
	elapsed = ((double)(end - start)) / CLOCKS_PER_SEC;
	fprintf(output_file, "\nTarget raggiunto, soluzione trovata in %f sec\n\n", elapsed);

	//print_matrix(matrix, row, col, output_file);
	print_incognite(matrix, row, col, numero_variabili, monomi, output_file);
	for (int i = 0; i<row*col; i++) {
		if (matrix[i] > module) {
			printf("OVERFLOW\n");
		}
	}

	free(matrix);
	free(degree);
	cudaDeviceReset();

	return 0;
}

