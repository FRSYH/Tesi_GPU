#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <time.h>
#include <omp.h>
#include <stdbool.h>
#include "linalg.h"
#include "matrix.h"
#include "scan.h"


//dichiarazione variabili globali
int max_degree = 0;
long long module = 0;
bool verbose_flag, help_flag, version_flag, test_flag, loop_flag, verify_flag, manual_expand_flag, reduced_expand_flag, set_expand_flag, partial_gauss_flag;

struct map_row {
	int len;
	int *col;
};

struct map {
	int len;
	struct map_row *row;
};


//moltiplica la riga indicata per tutti i possibili monomi
void moltiplica_riga(long long ***m, int *row, int col, int riga, struct map map,int * degree, int **vet, int num_var);

void moltiplica_riga_forn(long long ***m, int * row, int col, int riga, struct map map,int * degree, int **vet, int num_var, int stop_degree);

int init_matrix(long long **m, int row, int col,int **vet_grd, char *v, int num_var,int (*ord) (const void *, const void *, void*),FILE *input_file); //inizializza la matrice dei coefficienti

//initialize the vector that keeps the number of monomial with the same grade and their position
void init_degree_vector(int * degree, int num_var);


//return the grade of a monomial
int grado_monomio(int posizione, int **vet, int num_var);


//restituisce un array contenente tutti i len monomi con n variabili e grado <= m
//len è il numero di possibili monomi con n variabili e grado <= m
int **monomial_computation(int n, int m, int len);

//funzione ricorsiva che calcola tutti i possibili monomi con n variabili e grado <= m
//e li inserisce nell'array vet. Chiamata da monomial_computation
void monomial_computation_rec(int n, int m,  int **vet, int turn, int *monomial, int *pos);

//mappa tutte le possibili moltiplicazioni dei monomi di n variabili e grado <= m
//dell'array vet di lunghezza len, nella matrice map[len][len]. compar è la funzione
//secondo cui vet è ordinato.
void setup_map(int **map, int **vet, int len, int n, int m, int (*compar) (const void *, const void *, void*));

void setup_struct_map(struct map *map, int **vet, int len, int n, int m, int (*compar) (const void *, const void *, void*));

void print_struct_map(struct map map, FILE *output_file);

void free_struct_map(struct map *map);

//multiply the entire matrix for all possible correct monomial
void moltiplica_matrice(long long ***m, int *row, int col, struct map map, int *degree, int **vet, int num_var, int start);

//compute the different degree of the matrix rows
void matrix_degree(long long **m, int row, int col, int *m_deg, int **vet, int num_var);

//formatted print of matrix degree
void print_matrix_degree(int *m_deg, FILE *output_file);

//funzione di confronto gli array rowA con rowB, scorrendo gli elementi da destra a sinistra
//restituisce 1 se rowA > rowB, -1 se rowB > rowA, 0 altrimenti. Compatibile con qsort_r
int compare_arrays(const void *rowA, const void *rowB, void *columns);

//funzione che prende un vettore vet contenente length vettori di lunghezza (max_deg+1)
//retituisce il numero di cicli di fila partendo dal fondo di vet
int find_finishing_cycle(int **vet, int length, int max_deg);

//compare the matrix degree with the target degree wich are {0,1,2,3,4,5...max_degree} return 0 if equal -1 if not
int target_degree(int *v);

void partial_gauss(long long **m, int row, int col, int num_var, FILE *output_file);

int compare_arrays(const void *rowA, const void *rowB, void *columns);

void set_expansion_degree(int *expansion_degree,int *m_deg);

void execute_eliminazione(long long ***m, int * d_row, int col, struct map map, int *degree, int **vet, int num_var, int n_loops, int expansion, FILE *output_file);

void execute_confronto(long long ***m, int * d_row, int col, struct map map, int *degree, int **vet, int num_var, int n_loops, int expansion, FILE *output_file);

void execute_standard(long long ***m, int * d_row, int col, struct map map, int *degree, int **vet, int num_var, int n_loops, int expansion, FILE *output_file);

void verifica_correttezza(long long **m, int row, int col, struct map map, int *degree, int **vet, int num_var, int n_loops, int expansion, FILE *output_file, int ex1, int ex2);

void print_incognite(long long **m, int row, int col, int num_var, int **vet, FILE *output_file);


int main (int argc, char *argv[]){

	FILE *input_file = NULL, *output_file = NULL;
	//inizializzo flag a false
	verbose_flag = help_flag = version_flag = test_flag= loop_flag = verify_flag = manual_expand_flag = reduced_expand_flag = set_expand_flag = partial_gauss_flag = false;
	//parsing delle opzioni
	//numero cicli dei gradi di default = 30
	//execute_standard di default
	//verifica tra standard e confronto di default
	//0 -> standard
	//1 -> confronto
	//2 -> eliminazione
	int n_loops = 30, execute = 0, ex1 = 0, ex2 = 1, expansion = 1;

	for (int parsed = 1; parsed < argc; parsed++) {
		if (!strcmp(argv[parsed], "--verbose"))
			verbose_flag = true;
		else if (!strcmp(argv[parsed], "--help"))
			help_flag = true;
		else if (!strcmp(argv[parsed], "--test"))
			test_flag = true;
		else if (!strcmp(argv[parsed], "--manual-expand"))
			manual_expand_flag = true;
		else if (!strcmp(argv[parsed], "--reduced-expand"))
			reduced_expand_flag = true;
		else if (!strcmp(argv[parsed], "--partial-gauss"))
			partial_gauss_flag = true;
		else if (!strcmp(argv[parsed], "--loops")) {
			loop_flag = true;
			parsed++;
			if (atoi(argv[parsed]) > 0)
				n_loops = atoi(argv[parsed]);
		}
		else if (parsed < argc && !strcmp(argv[parsed], "--input")) {
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
		else if (parsed < argc && !strcmp(argv[parsed], "--execute")) {
			parsed++;
			execute = atoi(argv[parsed]);
		}
		else if (parsed < argc && !strcmp(argv[parsed], "--set-expand")) {
			set_expand_flag = true;
			parsed++;
			if (atoi(argv[parsed]) > 0)
				expansion = atoi(argv[parsed]);
		}
		else if (parsed < argc && !strcmp(argv[parsed], "--verify")) {
			verify_flag = true;
			parsed++;
			ex1 = atoi(argv[parsed]);
			parsed++;
			ex2 = atoi(argv[parsed]);
		}
	}

	if (help_flag) {
		fprintf(stdout, "Per istruzioni sull'utilizzo --> https://github.com/FRSYH/Tesi\n");
		return 0;
	}

	if (!input_file)
		input_file = stdin;
	if (!output_file)
		output_file = stdout;

	//inizio
	double start_time = omp_get_wtime(), stopwatch;

	int row, col, num_var,n;
	int *d_row,**vet, len, **map;
	long long **m;
	char *v;
	row = col = num_var = 0;
	int (*ord) (const void *, const void *, void*);

	struct map smap;


	//alloca la matrice principale, legge da input: il modulo,massimo grado e numero variabili
	allocation(&m,&row,&col,&num_var,&v,&n,&module,&max_degree,input_file);
	d_row = &row;

	if( order(&ord,n) != 0 ){
		fprintf(stderr, "Ordinamento insesistente!!!\n\nTERMINAZIONE PROGRAMMA");
		return 0;
	}


	int degree[max_degree+1];
	len = col;

	//crea il vettore con tutti i possibili monomi avendo num_var varaibili e max_degree come massimo grado
	vet = monomial_computation(num_var, max_degree, len);


	//ordina il vettore dei monomi secondo un determinato ordinamento, ordinamento intercambiabile
	qsort_r(vet, len, sizeof(int*), ord, &num_var);

	//inizializzazione matrice (lettura dati input)
	if( init_matrix(m,row,col,vet,v,num_var,ord,input_file) == -1 ){
		fprintf(stderr, "Errore di input !!!\n\nTERMINAZIONE PROGRAMMA"); //se l'input è in formato scorrettro abort del programma
		return 0;
	}

	fprintf(output_file, "\nInizializzazione in %f sec\n",omp_get_wtime()-start_time);
	stopwatch = omp_get_wtime();
	//allocazione matrice che mappa le posizioni dei prodotti dei monomi
	//matrix_alloc_int(&map,len,len);
	//creazione della mappa
	//setup_map(map, vet, len, num_var, max_degree,ord);

	setup_struct_map(&smap,vet, len, num_var, max_degree,ord);

	fprintf(output_file, "\nMappa creata in %f sec,   %d x %d \n\n",omp_get_wtime()-stopwatch,len,len);


	//RISOLUZIONE PROBLEMA
	
	//testing
	double t0 = omp_get_wtime();

	//inizializzazione vettore dei gradi dei polinomi
	init_degree_vector(degree,num_var);

	//eseguo moltiplicazione e riduzione di Gauss finche non trovo soluzione
//----------------------------------------------------------------------------

	if (!verify_flag)
		switch (execute) {
		case 1:
			execute_confronto(&m,d_row,col,smap,degree,vet,num_var,n_loops,expansion,output_file);
			break;
		case 2:
			execute_eliminazione(&m,d_row,col,smap,degree,vet,num_var,n_loops,expansion,output_file);
			break;
		default:
			execute_standard(&m,d_row,col,smap,degree,vet,num_var,n_loops,expansion,output_file);
			break;	
		}
	else	
		verifica_correttezza(m,row,col,smap,degree,vet,num_var,n_loops,expansion,output_file,ex1,ex2);

//----------------------------------------------------------------------------

	
	double t1 = omp_get_wtime();

	fprintf(output_file, "\nTarget raggiunto, soluzione trovata in %f sec\n\n",omp_get_wtime()-start_time);


	//stampa tutta la matrice soluzione
	if (verbose_flag) {
		fprintf(output_file, "\n\n matrice soluzione:\n\n");
		print_matrix(m, row, col, output_file);
	}
	fprintf(output_file, "Valori delle incognite\n\n");
	print_incognite(m,row,col,num_var,vet,output_file);

	//deallocazione di tutti i puntatori utilizzati
	matrix_free_long(&m,row,col);
	free_struct_map(&smap);
	//matrix_free_int(&map,len,len);
	matrix_free_int(&vet,len,num_var);

	//testing
	if (test_flag) {fprintf(output_file, "\n%f\n", t1-t0);}

	return 0;
}





void execute_eliminazione(long long ***m, int * d_row, int col, struct map map, int *degree, int **vet, int num_var, int n_loops, int expansion, FILE *output_file){

	//eseguo moltiplicazione e riduzione di Gauss finche non trovo soluzione
	//non moltiplico le linee iniziali uguali a quelle dell'iterazione precedente
	//tot = matrice totale su cui faccio gauss
	//prev = matrice dell'iterazione precedente
	//m = "now" = matrice che contiene le righe diverse tra tot e prev, che vanno moltiplicate
	//	e aggiunte a tot prima di fare gauss
	double start_time = omp_get_wtime(), stopwatch;

	//creo l'array che conterrà i gradi dei vari round
	int **m_deg_array = malloc(sizeof(int*));
	m_deg_array[0] = calloc(max_degree+1, sizeof(int));
	int n_round = 0;
	int *m_deg = m_deg_array[0];
	
	fprintf(output_file, "Inizio computazione, metodo eliminazione\n");
	matrix_degree(*m,*d_row,col,m_deg,vet,num_var);

	int flag,old,new;
	flag = old = new = 0;
	old = *d_row;
	//assumo espansione normale
	int expansion_degree = max_degree;

	long long **prev = NULL;
	int row_prev = 0, row_tot = *d_row;
	long long **tot = NULL;
	//tot = now
	matrix_alloc_long(&tot, row_tot, col);
	matrix_cpy(*m, row_tot, col, tot);

	if (set_expand_flag)
		expansion_degree = expansion;


	while( flag != 1 ){
		n_round++;
		
		//prev = tot, salvo la matrice della precedente iterazione
		row_prev = row_tot;
		matrix_alloc_long(&prev, row_prev, col);
		matrix_cpy(tot, row_prev, col, prev);

		//mult(now), moltiplico le linee diverse. Nella prima iterazione
		//non sono diverse, ma sono poche e non influisce sulle prestazioni
		fprintf(output_file, "\n -Eseguo moltiplicazione su m, ");
		fflush(stdout);
	
		set_expansion_degree(&expansion_degree, m_deg);

		int length = *d_row;
		stopwatch = omp_get_wtime();
		for(int i=0; i<length; i++){
			moltiplica_riga_forn(m,d_row,col,i,map,degree,vet,num_var,expansion_degree);	
		}

		//moltiplica_matrice(m,d_row,col,map,degree,vet,num_var,0);
		
		//tot = tot + now, faccio append a tot di now (linee moltiplicate)
		append_and_free_matrix(&tot, &row_tot, col, m, d_row, col);
		//non mi serve più now
		//matrix_free_long(m, *d_row, col);
		fprintf(output_file, "grado di espansione: %d, numero righe: %d     (%f sec)\n", expansion_degree, row_tot, omp_get_wtime()-stopwatch);

		if (partial_gauss_flag)
			partial_gauss(tot, row_tot, col, num_var, output_file);
		
		//gauss(tot), gauss su tot che ha anche le linee moltiplicate
		fprintf(output_file, " -Eseguo Gauss, ");
		fflush(stdout);
		stopwatch = omp_get_wtime();	
		
		gauss(tot, row_tot, col, module, 0, NULL);
		eliminate_null_rows(&tot, &row_tot, col);
		fprintf(output_file, "numero righe: %d              (%f sec)\n", row_tot ,omp_get_wtime()-stopwatch);
		
		//now = tot - prev, tolgo da tot le linee iniziali uguali a quelle
		//dell'iterazione precedente e assegno il risultato a now
		*d_row = row_tot;
		matrix_alloc_long(m, *d_row, col);
		matrix_cpy(tot, row_tot, col, *m);
		//eliminate_equal_rows(m, d_row, prev, row_prev, col);
		eliminate_equal_starting_rows(m, d_row, prev, row_prev, col);
		
		//aggiungo all'array i gradi dell'attuale round
		//n_round+1 perchè salvo anche i gradi prima di inziare i round
		m_deg_array = realloc(m_deg_array, sizeof(int*)*(n_round+1));
		m_deg_array[n_round] = calloc(max_degree+1, sizeof(int));
		m_deg = m_deg_array[n_round];
		
		//degree(tot), aggiorno gradi/target
  		matrix_degree(tot, row_tot,col,m_deg,vet,num_var);
		print_matrix_degree(m_deg,output_file);
		new = row_tot;
		
		matrix_free_long(&prev, row_prev, col);


		//se i gradi dei round formano più di n_loops cilcli e se il flag è true
		//o trovo una matrice con le stesso numero di righe della precedente mi fermo
		if( (find_finishing_cycle(m_deg_array, n_round+1, max_degree) > n_loops && loop_flag) || ((!loop_flag) && old == new)  ) {
			flag = 1;
			fprintf(output_file, "\n\nEXIT: superato numero di cicli massimo o numero righe rimaste invariate\n\n");	
		}
		else
			if( target_degree(m_deg) == 0 )
				flag = 1;
			else{
				old = new;
				//verbose
				if (verbose_flag) {
					fprintf(output_file, "\nMatrice intermedia:\n\n");
					print_matrix(tot, row_tot, col, output_file);
				}
			}

	}

	matrix_free_long(m, *d_row, col);
	for (int i = 0; i < n_round+1; i++)
		free(m_deg_array[i]);	
	free(m_deg_array);
	*m = tot;
	*d_row = row_tot;
}



void execute_confronto(long long ***m, int * d_row, int col, struct map map, int *degree, int **vet, int num_var, int n_loops, int expansion, FILE *output_file){

	
	int flag,old,new,inv;
	flag = old = new = 0;
	old = *d_row;

	int st = inv = 0;
	
	int *v1,*v2;

	double start_time = omp_get_wtime(), stopwatch;

	//creo l'array che conterrà i gradi dei vari round
	int **m_deg_array = malloc(sizeof(int*));
	m_deg_array[0] = calloc(max_degree+1, sizeof(int));
	int n_round = 0;
	int *m_deg = m_deg_array[0];
	//assumo espansione normale
	int expansion_degree = max_degree;
	matrix_degree(*m,*d_row,col,m_deg,vet,num_var);

	if (set_expand_flag)
		expansion_degree = expansion;

	fprintf(output_file, "Inizio computazione, metodo confronto\n");
	//-------------------------------------------------------------------------------------------

	fprintf(output_file, "\n -Eseguo moltiplicazione, ");
	fflush(stdout);

	stopwatch = omp_get_wtime();

	set_expansion_degree(&expansion_degree, m_deg);

	//moltiplico
	int length = *d_row;
	for(int i=0; i<length; i++){
		moltiplica_riga_forn(m,d_row,col,i,map,degree,vet,num_var,expansion_degree);	
	}


	//moltiplico la matrice per tutti i monomi possibili
	//moltiplica_matrice(m,d_row,col,map,degree,vet,num_var,0);
	
	fprintf(output_file, "grado di espansione: %d, numero righe: %d     (%f sec)\n", expansion_degree, *d_row, omp_get_wtime()-stopwatch);

	while( flag != 1 ){
		n_round++;
//-------------------------------------------------------------------------------------------
		// calcolo la posizone dell'ultimo elemento di ogni riga della matrice prima di effettuare gauss	

		v1 = calloc( *d_row , sizeof( int ) );
		for( int i = 0; i < *d_row; i++ ){
			for( int j = col-1; j>=0; j--){
				if( (*m)[i][j] != 0 ){
					v1[i] = j;
					break;	
				}	
			}
		}

		//passo il vettore appena calcolato alla procedura di gauss per invertire le righe in modo analogo a quanto avviene nella riduzione

//-------------------------------------------------------------------------------------------
		if (partial_gauss_flag)
			partial_gauss(*m, *d_row, col, num_var, output_file);


		fprintf(output_file, "\n -Eseguo Gauss, ");
		fflush(stdout);
		stopwatch = omp_get_wtime();

		//applico la riduzione di Gauss
		gauss(*m, *d_row, col, module, st,v1);
		//magma_gauss(m, *d_row, col, module);

		//elimino le righe nulle della matrice
		eliminate_null_rows(m,d_row,col);

		//aggiungo all'array i gradi dell'attuale round
		//n_round+1 perchè salvo anche i gradi prima di inziare i round
		m_deg_array = realloc(m_deg_array, sizeof(int*)*(n_round+1));
		m_deg_array[n_round] = calloc(max_degree+1, sizeof(int));
		m_deg = m_deg_array[n_round];

		fprintf(output_file, "numero righe: %d               (%f sec)\n", *d_row,omp_get_wtime()-stopwatch);
  		matrix_degree(*m,*d_row,col,m_deg,vet,num_var);
		print_matrix_degree(m_deg,output_file);

		new = *d_row;
		st = new;
//-----------------------------------------------------------
		// ricalcolo la posizione dell'ultimo elemento di ogni riga dopo aver effettuato gauss
			
		v2 = calloc( *d_row , sizeof( int ) );
		for( int i = 0; i < *d_row; i++ ){
			for( int j = col-1; j>=0; j--){
				if( (*m)[i][j] != 0 ){
					v2[i] = j;
					break;	
				}	
			}
		}
//----------------------------------------------------------
		// controllo le condizioni di uscita

		//se i gradi dei round formano più di n_loops cilcli e se il flag è true
		//o trovo una matrice con le stesso numero di righe della precedente mi fermo
		if( (find_finishing_cycle(m_deg_array, n_round+1, max_degree) > n_loops && loop_flag) || (!(loop_flag) && old == new)  ) {
			flag = 1;
			fprintf(output_file, "\n\nEXIT: superato numero di cicli massimo o numero righe rimaste invariate\n\n");	
			break;
		}else{
			if( target_degree(m_deg) == 0 ){
				flag = 1;
				break;
			}
			else{
				old = new;
				//verbose
				if (verbose_flag) {
					fprintf(output_file, "\nMatrice intermedia:\n\n");
					print_matrix(*m, *d_row, col, output_file);
				}
			}
		}
		fprintf(output_file, "\n -Eseguo moltiplicazione, ");
		fflush(stdout);

		set_expansion_degree(&expansion_degree, m_deg);
		
		//moltiplico
		int length = *d_row;
		stopwatch = omp_get_wtime();
		for(int i=0; i<length; i++){
			if( v1[i] > v2[i] ){ 	//significa che la riga è stata ridotta

				moltiplica_riga_forn(m,d_row,col,i,map,degree,vet,num_var,expansion_degree);	//allora moltiplico tale riga
				//moltiplica_riga(m,d_row,col,i,map,degree,vet,num_var);	//allora moltiplico tale riga
			}
		}

		fprintf(output_file, "grado di espansione: %d, numero righe: %d     (%f sec)\n", expansion_degree, *d_row, omp_get_wtime()-stopwatch);

	
		free(v2);
		free(v1);
	}
	for (int i = 0; i < n_round+1; i++)
		free(m_deg_array[i]);	
	free(m_deg_array);
	//finito algoritmo moltiplicazione e riduzione
}




void execute_standard(long long ***m, int * d_row, int col, struct map map, int *degree, int **vet, int num_var, int n_loops, int expansion, FILE *output_file){

	double start_time = omp_get_wtime(), stopwatch;

	//creo l'array che conterrà i gradi dei vari round
	int **m_deg_array = malloc(sizeof(int*));
	m_deg_array[0] = calloc(max_degree+1, sizeof(int));
	int n_round = 0;
	int *m_deg = m_deg_array[0];

	fprintf(output_file, "Inizio computazione, metodo standard\n");
	matrix_degree(*m,*d_row,col,m_deg,vet,num_var);

	int flag,old,new;
	flag = old = new = 0;
	old = *d_row;
	//assumo espansione normale
	int expansion_degree = max_degree;
	int st = 0;

	
	if (set_expand_flag)
		expansion_degree = expansion;


	while( flag != 1 ){
		n_round++;

		fprintf(output_file, "\n -Eseguo moltiplicazione, ");
		fflush(stdout);

		set_expansion_degree(&expansion_degree, m_deg);
		stopwatch = omp_get_wtime();
		
		//moltiplico
		int length = *d_row;
		for(int i=0; i<length; i++){
			moltiplica_riga_forn(m,d_row,col,i,map,degree,vet,num_var,expansion_degree);	
		}

		//moltiplico la matrice per tutti i monomi possibili
		//moltiplica_matrice(m,d_row,col,map,degree,vet,num_var,0);
		
		fprintf(output_file, "grado di espansione: %d, numero righe: %d     (%f sec)\n", expansion_degree, *d_row, omp_get_wtime()-stopwatch);

		if (partial_gauss_flag)
			partial_gauss(*m, *d_row, col, num_var, output_file);

		fprintf(output_file, "\n -Eseguo Gauss, ");
		fflush(stdout);
		stopwatch = omp_get_wtime();	


		//applico la riduzione di Gauss
		gauss(*m, *d_row, col, module, st,NULL);
		//elimino le righe nulle della matrice
		eliminate_null_rows(m,d_row,col);
		
		
		//aggiungo all'array i gradi dell'attuale round
		//n_round+1 perchè salvo anche i gradi prima di inziare i round
		m_deg_array = realloc(m_deg_array, sizeof(int*)*(n_round+1));
		m_deg_array[n_round] = calloc(max_degree+1, sizeof(int));
		m_deg = m_deg_array[n_round];
		
		fprintf(output_file, "numero righe: %d               (%f sec)\n", *d_row,omp_get_wtime()-stopwatch);
  		matrix_degree(*m,*d_row,col,m_deg,vet,num_var);
		print_matrix_degree(m_deg,output_file);

		new = *d_row;
		st = new;


		//se i gradi dei round formano più di n_loops cilcli e se il flag è true
		//o trovo una matrice con le stesso numero di righe della precedente mi fermo
		if( (find_finishing_cycle(m_deg_array, n_round+1, max_degree) > n_loops && loop_flag) || (!(loop_flag) && old == new)  ) {
			flag = 1;
			fprintf(output_file, "\n\nEXIT: superato numero di cicli massimo o numero righe rimaste invariate\n\n");	
		}
		else
			if( target_degree(m_deg) == 0 )
				flag = 1;
			else{
				old = new;
				//verbose
				if (verbose_flag) {
					fprintf(output_file, "\nMatrice intermedia:\n\n");
					print_matrix(*m, *d_row, col, output_file);
				}
			}

	}
	for (int i = 0; i < n_round+1; i++)
		free(m_deg_array[i]);	
	free(m_deg_array);

//finito algoritmo moltiplicazione e riduzione
}


int init_matrix(long long **m, int row, int col, int **vet_grd, char *v, int num_var, int (*ord) (const void *, const void *, void*), FILE *input_file ){
//Inizializza la matrice principale (dei coefficienti) con i coefficienti dei polinomi forniti come input.
	return parse(num_var,v,m,row,vet_grd,col,module,ord,input_file);
}


/*
Moltiplica la riga indicata (riga) della matrice m per ogni monomio con grado <= al stop_degree
o comunque il polinomio risultante abbia grado <= max_degree.
Il prodotto avviene su monomi con coefficiente sempre uguale a 1.
Il prodotto consiste quindi in uno shift della posizione del monomio in esame.
Il parametro map fornisce una mappa delle posizioni in cui inserire il prodotto di due monomi.
La matrice aumenta il numero di righe in base a quanti prodotti devo eseguire.
Gli ultimi due parametri servono per il calcolo del grado di un monomio.
*/

void moltiplica_riga_forn(long long ***m, int * row, int col, int riga, struct map map, int * degree, int **vet, int num_var, int stop_degree){

	int grado_massimo_riga, grado_massimo_monomio,i,j,last,new_row;
	last = -1;
	//cerco la posizione dell'ultimo coefficiente non nullo del polinomio rappresentato nella riga.
	for(i=col-1; i>0; i--){
		if( (*m)[riga][i] != 0 ){
			last = i;
			break;
		}
	}
	//risalgo al grado del monomio appena trovato
	//scorro la lista delle posizioni di inizio dei monomi con lo stesso grado
	if( last != -1 ){

		grado_massimo_riga = grado_monomio(last,vet,num_var);

		//calcolo il grado massimo che deve avere il monomio per cui moltiplicare
		grado_massimo_monomio = max_degree - grado_massimo_riga;
		// a questo punto conosco per quanti monomi devo moltiplicare e quindi
		// conosco il numero di righe che devo aggiungere alla matrice
		new_row = 0;

		if( stop_degree != 0 ){
			
			if( grado_massimo_monomio > stop_degree ){
				grado_massimo_monomio = stop_degree;
			}
		}

		for(i=1; i<(grado_massimo_monomio+1); i++){
			new_row += degree[i];
		}

		*m = realloc( *m , (*row + new_row ) * sizeof (long long *));

		for (i=(*row); i< (*row + new_row ); i++)
			(*m)[i] = calloc(col , sizeof (long long) );


		for(i=1; i<(new_row+1); i++){     								//scorre tutti i monomi per i quali posso moltiplicare
			//#pragma omp parallel for shared (m,row,riga)
			for(j=0; j<(last+1); j++){     								//scorre fino all'ultimo elemento della riga
				(*m)[*row][ map.row[i].col[j] ] = (*m)[riga][j];  				//shift nella posizione corretta indicata dalla mappa
			}
			*row = *row + 1;											//aumento del conteggio delle righe
		}
	}
}




int grado_monomio(int posizione, int **vet, int num_var){
//Calcola il grado del monomio a partire dalla posizione occupata nel vettore (ordinato) delle posizioni rispetto l'ordinamento scelto.
//(la posizione occupata deve essere corretta).
	int i,grado;
	grado = 0;
	for(i=0; i<num_var; i++){
		grado += vet[posizione][i];
	}
	return grado;
}


void init_degree_vector(int * degree, int num_var){
//inizializza il vettore degree con il numero di monomi di grado i-esimo <= del grado massimo
	int i,j,c;
	for(i=0; i<max_degree+1; i++){
		c = combination(num_var,i);
		degree[i] = c;
	}
}



/*
mappa tutte le possibili moltiplicazioni dei monomi di n variabili e grado <= m
dell'array vet di lunghezza len, nella matrice map[len][len].
Al termine map[x][y] contiene la posizione all'interno di vet del
monomio risultato dal prodotto di vet[x]*vet[y]
Esempio: vet[4] * vet[10] = vet [map[4][10]]
se il grado del prodotto supera m viene messo il valore -1 nella matrice.
compar è la funzione secondo cui vet è ordinato,
la matrice map deve essere già correttamente allocata
*/
void setup_map(int **map, int **vet, int len, int n, int m, int (*compar) (const void *, const void *, void*)) {

	int sum, *temp = malloc(n * sizeof(int));

	//per ogni monomio in vet

	for (int row = 0; row < len; row++)
		//provo a moltiplicarlo con ogni monomio in vet
		for (int col = 0; col < len; col++) {
			sum = 0;
			//eseguo il prodotto (sum è la somma dei gradi)
			for (int v = 0; v < n; v++) {
				temp[v] = vet[row][v] + vet[col][v];
				sum += temp[v];
			}
			//se il grado del prodotto > grado massimo tutti i restanti prodotti
			//su quella riga sono > grado massimo, setto a -1 il resto della riga
			if (sum > m) {

				for (int i = col; i < len; i++)
					map[row][i] = -1;
				break;
			}
			//altrimenti cerco il prodotto in vet e metto l'indice in map
			else{
				map[row][col] = (int **)(bsearch_r((void *) &temp, (void *) vet, len, (sizeof(int*)), compar, &n)) - vet;
			}
		}
	free(temp);
}


void setup_struct_map(struct map *map, int **vet, int len, int n, int m, int (*compar) (const void *, const void *, void*)  ){
	
	int sum, index=len;

	//	inizializzo la struttura map, la mappa ha len righe.
	map->len = len;
	map->row = malloc( map->len * sizeof(struct map_row) );

	//per ogni monomio in vet
	int row, col, i, v;
	#pragma omp parallel for private(row,col,sum,index,i,v) schedule(dynamic)
	for (row = 0; row < len; row++){
		index = 0;
		//dichiarati dentro per la parallelizzazione
		int *temp = malloc(n * sizeof(int));
		int *save = calloc(len, sizeof(int));
		//provo a moltiplicarlo con ogni monomio in vet
		for (col = 0; col < len; col++) {
			sum = 0;
			//eseguo il prodotto (sum è la somma dei gradi)
			for (v = 0; v < n; v++) {
				temp[v] = vet[row][v] + vet[col][v];
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
			else{
				save[col] = (int **)(bsearch_r((void *) &temp, (void *) vet, len, (sizeof(int*)), compar, &n)) - vet;
			}
		}

		//	terminato il ciclo sulle colonne posso inizializzare la struttura perchè conosco tutti gli elementi da inserire	
		//  la riga attuale ha esattamente index elementi diversi da -1, quindi la riga avrà lunghezza pari a index precedentemente calcolato
		//  alloco la riga con un array da index elementi

		map->row[row].len = index;
		map->row[row].col = malloc( map->row[row].len * sizeof(int) );
		//	a questo map devo copiare gli elementi generati dento alla struttura

		for(i=0; i<map->row[row].len; i++)
			map->row[row].col[i] = save[i];
		
		free(temp);
		free(save);
	}
}

void print_struct_map(struct map map, FILE *output_file){

	fprintf(output_file, "Inizio stampa\n");
	for(int i=0; i<map.len; i++ ){
		for(int j=0; j<map.row[i].len; j++){
			fprintf(output_file, "%d ", map.row[i].col[j]);
		}
		fprintf(output_file, "\n");
	}
	fprintf(output_file, "Fine stampa\n");
}


void free_struct_map(struct map *map){
	for(int i=0; i<map->len; i++){
		free(map->row[i].col);
	}
	free(map->row);
}


/*restituisce un array contenente tutti i len monomi con n variabili e grado <= m
len è il numero di possibili monomi con n variabili e grado <= m
i monomi sono array di interi di lunghezza n dove il valore di ogni posizione rappresenta
il grado della variabile in quella posizione. Esempio: n=3, x^2*y*z = [2,1,1]
len viene passato come argomento per evitare di ricalcolarlo internamente
*/
int **monomial_computation(int n, int m, int len) {

	int **vet, *monomial;

	//alloco la memoria per l'array
	matrix_alloc_int(&vet,len,n);

	//strutture di supporto necessarie per il calcolo
	monomial = malloc(n * sizeof(int));
	int pos = 0;

	//il calcolo è fatto dalla funzione ricorsiva correttemente parametrizzata
	monomial_computation_rec(n, m, vet, 0, monomial, &pos);

	free(monomial);

	return vet;
}



/*funzione ricorsiva che calcola tutti i possibili monomi con n variabili e grado <= m
e li inserisce nell'array vet. I monomi sono rappresentati come array di interi dove
il valore di ogni posizione rappresenta il grado della variabile in quella posizione.
Esempio: n=3, x^2*y*z = [2,1,1].
L'array vet deve essere già allocato correttamente. Gli altri parametri sono necessari
per la struttura ricorsiva della funzione e alla prima chiamata devono essere:
- turn = 0, rappresenta la posizione della variabile nel monomio
- monomial = array di interi di lunghezza n già allocato e usato per calcolare i vari monomi
- *pos = 0 puntatore ad intero, rappresenta la prima posizione libera nell'array vet
*/
void monomial_computation_rec(int n, int m, int **vet, int turn, int *monomial, int *pos) {

	//per ogni variabile provo tutti i gradi da 0 a m
	for (int degree = 0; degree <= m; degree++) {
		//se questa è la prima variabile azzero il monomio
		if (turn == 0) {
			//azzero il monomio lasciando solo il grado della prima variabile
			monomial[0] = degree;
			for (int v = 1; v < n; v++)
				monomial[v] = 0;
		}
		//altrimenti le altre variabili aggiungo il proprio grado al monomio
		else
			monomial[turn] = degree;


		//ottengo il grado del monomio sommando i gradi delle variabili
		int sum = 0;
		for (int v = 0; v <= turn; v++)
			sum += monomial[v];
		//se il grado del monomio supera quello massimo non ha senso continuare a cercare
		//altri monomi partendo da questo, perchè tutti avranno grado maggiore o uguale
		if (sum > m)
			break;

		//se questa è l'ultima variabile copia il monomio nell'array vet and incrementa l'indice pos
		if (turn == (n-1)) {
			vctcpy(vet[(*pos)], monomial, n);
			(*pos)++;
		}
		//altrimenti richiama se stessa cambiando la variabile (turn)
		else
			monomial_computation_rec(n, m, vet, turn+1, monomial, pos);
	}

	return;
}




void matrix_degree(long long **m, int row, int col, int *m_deg, int **vet, int num_var){
//m_deg è un vettore che ha lunghezza pari al grado massimo.
//la funzione calcola i gradi dei polinomi presenti nella matrice.
//Ogni cella del vettore m_deg rappresenta un grado, se esso compare nella matrice allora viene impostato a 1 o altrimenti.

	int i,j,last,grado;
	for(i=0; i<row; i++){
		for(j=col-1; j>0; j--){
			if( m[i][j] != 0 ){
				last = j;           //posizione dell'ultimo coefficiente della riga
				break;
			}
		}
		grado = grado_monomio(last,vet,num_var);
		m_deg[grado] = 1;
	}
}

//in base alle flag calcola il grado di espansione e lo metto in expansion_degree
//m_deg è l'array dei gradi dei polinomi presenti nella matrice
void set_expansion_degree(int *expansion_degree,int *m_deg) {

	if (manual_expand_flag) {
		fprintf(stdout, "inserire grado di espansione: ");
		fscanf(stdin, " %d", expansion_degree);
	}
	//se ridotta calcolo il grado mancante
	else if (reduced_expand_flag)
		for(int i=max_degree; i>0; i--){
			if( m_deg[i] == 0 ){
				*expansion_degree = i;
				break;
			}
		}
	return;

}


//esegue gauss, se richiesto, su una porzione della matrice specificata dall'utente
//e stampa il risultato nel output_file
void partial_gauss(long long **m, int row, int col, int num_var, FILE *output_file) {

	char answer = 'n';
	long long **temp;
	int col_degree, temp_col;

	fprintf(stdout, "\nSi vuole effettuare gauss su una porzione della matrice? (y/n):  ");
	fscanf(stdin, " %c", &answer);
	
	if (answer != 'y')
		return;

	fprintf(stdout, "\nScegliere il grado massimo dei monomi su cui eseguire l'eliminazione:  ");
	fscanf(stdin, " %d", &col_degree);

	if (col_degree > max_degree || col_degree < 0)
		col_degree = max_degree;
	temp_col = monomial_combinations(num_var, col_degree);

	matrix_alloc_long(&temp, row, temp_col);
	matrix_cpy(m, row, temp_col, temp);
	gauss(temp, row, temp_col, module, 0, NULL);
	eliminate_null_rows(&temp, &row, temp_col);
	
	fprintf(output_file, "\nEseguo gauss parziale su %d colonne\n\n", temp_col);
	print_matrix(temp, row, temp_col, output_file);

	matrix_free_long(&temp, row, temp_col);

	return;
}


void print_matrix_degree(int *m_deg, FILE *output_file){
//stampa il vettore dei gradi della matrice.
	int i;
	fprintf(output_file, "Gradi della matrice = {");
	for(i=0; i<max_degree+1; i++)
		if( m_deg[i] != 0 )	fprintf(output_file, " %d ",i);
	fprintf(output_file, "}\n");
}




int target_degree(int *v){
//Controlla se il vettore v rappresenta la condizione di terminazione con gradi completi {1,2,3,...,max_degree}
//Se la condizione è soddisfatta return 0 altrimenti -1

	int i,flag;
	flag = 0;
	for(i=1; i<max_degree+1; i++){
		if( v[i] != 1 ){
			flag = -1;
			break;
		}
	}
	return flag;
}


//funzione che prende un vettore vet contenente length vettori di lunghezza (max_deg+1)
//retituisce il numero di cicli massimo all'interno di vet partendo dai valori in fondo
int find_finishing_cycle(int **vet, int length, int max_deg) {
	int st1, st2, temp1, temp2, max_cycles, n_cycles;
	st1 = length - 1;
	max_cycles = n_cycles = 0;
	int l_deg = max_deg + 1;
	
	//per ogni elemento prima di st1 fino a metà
	//oltre la metà non ci può più essere un loop
	for (st2 = length-2; st2 >= (length/2)-1; st2--) {
		n_cycles = 0;
		//se sono uguali potrebbe essere l'inizio di un ciclo
		if (!compare_arrays(&(vet[st1]), &(vet[st2]), &l_deg)) {
			//guardo i valori precedenti agli start
			temp1 = st1 - 1;
			temp2 = st2 -1;
			//condizione necessaria per i primi valori di fila uguali
			//(fatto solo la prima volta) 12333 <- trovo i 333 = 2 cicli
			if (st1-1 == st2) {
				n_cycles++;
				while(temp2 >= 0 && !compare_arrays(&vet[st1], &vet[temp2], &l_deg)) {
					n_cycles++;
					temp2--;
				}
			}
			else
				while (temp1 > st2 && temp2 >= 0) {
					//se sono diversi non è un ciclo, esco
					if (compare_arrays(&vet[temp1], &vet[temp2], &l_deg))
						break;
					temp1--;
					temp2--;
					//sono arrivato a st2, ho fatto un ciclo e incremento
					if (temp1 == st2) {
						n_cycles++;
						//se non sono arrivato alla fine
						if (temp2 > 0) {
							st2 = temp2;
							temp1 = st1-1;
							temp2 = st2-1;
						}
					}
				}
				//aggiorno il numero di cicli massimo
				if (n_cycles > max_cycles)
					max_cycles = n_cycles;
		}
	}
	
	return max_cycles;
}

//funzione di confronto gli array rowA con rowB, scorrendo gli elementi da destra a sinistra
//restituisce 1 se rowA > rowB, -1 se rowB > rowA, 0 altrimenti. Compatibile con qsort_r
int compare_arrays(const void *rowA, const void *rowB, void *columns) {

	int *row1, *row2;
	int col;
	
	col = *((int *) columns);
	row1 = *((int **) rowA);
	row2 = *((int **) rowB);
	
	for (int i = col-1; i >= 0; i--) {
		if (row1[i] > row2[i])
			return 1;
		else if (row1[i] < row2[i])
			return -1;
	}
	
	return 0;
}

void verifica_correttezza(long long **m, int row, int col, struct map map, int *degree, int **vet, int num_var,int n_loops, int expansion, FILE *output_file, int ex1, int ex2){

	long long **m1,**m2,**m3;
	int row1,row2,row3;
	double start_time = omp_get_wtime(), stopwatch;

	row1 = row2 = row;

	matrix_alloc_long(&m1,row1,col);
	matrix_alloc_long(&m2,row2,col);

	matrix_cpy(m,row,col,m1);
	matrix_cpy(m,row,col,m2);



	fprintf(output_file, "\nESEGUO PRIMA RISOLUZIONE\n");
	switch (ex1) {
	case 1:
		execute_confronto(&m1,&row1,col,map,degree,vet,num_var,n_loops,expansion,output_file);
		break;
	case 2:
		execute_eliminazione(&m1,&row1,col,map,degree,vet,num_var,n_loops,expansion,output_file);
		break;
	default:
		execute_standard(&m1,&row1,col,map,degree,vet,num_var,n_loops,expansion,output_file);
		break;	
	}

	fprintf(output_file, "\nTERMINATA PRIMA RISOLUZIONE, NUMERO RIGHE:%d\n",row1);


	fprintf(output_file, "\n\nESEGUO SECONDA RISOLUZIONE\n");

	switch (ex2) {
	case 1:
		execute_confronto(&m2,&row2,col,map,degree,vet,num_var,n_loops,expansion,output_file);
		break;
	case 2:
		execute_eliminazione(&m2,&row2,col,map,degree,vet,num_var,n_loops,expansion,output_file);
		break;
	default:
		execute_standard(&m2,&row2,col,map,degree,vet,num_var,n_loops,expansion,output_file);
		break;	
	}


	fprintf(output_file, "\nTERMINATA SECONDA RISOLUZIONE, NUMERO RIGHE:%d\n",row2);

	fprintf(output_file, "\n\nCONFRONTO LE MATRICI PER OTTENERE IL NUMERO DI LINEE DIVERSE\n");

	//trovo il numero di linee diverse tra m1 e m2
	row3 = row1;
	matrix_alloc_long(&m3, row3, col);
	matrix_cpy(m1, row3, col, m3);
	eliminate_equal_rows(&m3, &row3, m2, row2, col);
	
	fprintf(output_file, "\nRIGHE DIVERSE RISPETTO TRA PRIMO E SECONDO METODO:%d\n", row3);
	matrix_free_long(&m3, row3, col);

	//controllo che i polinomi di una matrice siano linearmente
	//dipendenti da quelli dell'altra, necessario per soluzioni equivalenti
	append_and_free_matrix(&m1, &row1, col, &m2, &row2, col);

	fprintf(output_file, "\n\nESEGUO APPEND, NUMERO RIGHE:%d\n",row1);

	int *m_deg = calloc(max_degree+1, sizeof(int));
	matrix_degree(m1,row1,col,m_deg,vet,num_var);
	print_matrix_degree(m_deg,output_file);

	fprintf(output_file, " -Eseguo Gauss su matrice totale, ");
	fflush(stdout);
	stopwatch = omp_get_wtime();	
	
	gauss(m1, row1, col, module, 0, NULL);
	eliminate_null_rows(&m1, &row1, col);
	fprintf(output_file, "numero righe: %d              (%f sec)\n", row1 ,omp_get_wtime()-stopwatch);
	matrix_degree(m1,row1,col,m_deg,vet,num_var);
	print_matrix_degree(m_deg,output_file);	

	matrix_free_long(&m1,row1,col);

}


void print_incognite(long long **m, int row, int col, int num_var, int **vet, FILE *output_file){

	int grado,last;

	for(int r = row - (num_var+1); r<row; r++){

		//cerca la posizione dell'ulitmo elemento non nullo della riga r
		for( int i=col-1; i>=0; i-- ){
			if( m[r][i] != 0 ){
				last = i;
				break;
			}
		}
		//calcola il grado della riga r
		grado = grado_monomio(last,vet,num_var);
		//se il grado della riga r è 1 allora stampa tutta la riga della matrice
		if( grado == 1 ){
			for( int j=0; j<last+1; j++ ){
				fprintf(output_file, "%lli ", m[r][j]);
			}
			fprintf(output_file, "\n\n");
		}
	}
	fprintf(output_file, "\n");	
}
