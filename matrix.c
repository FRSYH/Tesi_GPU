#include <stdlib.h>
#include <stdio.h>
#include "matrix.h"

//Scambio di due righe della matrice m.
void swap_rows(long long **m, int row, int col, int j, int i){
	
	int k;
	long long tmp;
	if( j!=i ){
		#pragma omp parallel for private (k,tmp) shared (m)
		for(k=0;k<col;k++){
			tmp = m[i][k];
			m[i][k] = m[j][k];
			m[j][k] = tmp;
		}
	}
}

//Stampa formattata della matrice
void print_matrix(long long **m, int row, int col, FILE *output_file){

	int i,j;	
	for (i=0;i<row;i++)
	{
		for (j=0;j<col;j++){
			fprintf(output_file, "%lli ",m[i][j]);
		}
		fprintf(output_file, "\n\n\n");
	}
	fprintf(output_file, "\n");
}

//copia il vettore vet2 in vet1, entrambi di lunghezza len
void vctcpy(int *vet1, const int *vet2, int len) {
	for (int i = 0; i < len; i++)
		vet1[i] = vet2[i];
	return;
}

int null_rows(long long **m, int row, int col){
//calcola il numero di righe nulle presenti nella matrice m.

	int i,j,last,null_rows;
	null_rows = 0;
	for(i=0; i<row; i++){
		last = -1;
		for(j=col-1; j>-1; j--){
			if(m[i][j] != 0 ){
				last = j;
				break;
			}
		}
		if( last == -1 )
			null_rows++;
	}
	return null_rows;
}

void eliminate_null_rows(long long ***m, int *row, int col){
//Elimina dalla matrice m le righe nulle.
//N.B. questa procedura elimina le ultime righe nulle della matrice.
//Questa funzione DEVE essere utilizzata dopo la riduzione di Gauss.
//La riduzione di Gauss sposta nelle ultime posizioni tutte le righe nulle.
//Se non si esegue questa funzione dopo Gauss si possono eliminare righe non nulle.	

	int null_row = null_rows(*m,*row,col);
	int new_rows = *row - null_row;
	if(null_row != 0){
		for(int i=new_rows; i < *row; i++){
			free( (*m)[i] );
		}
		*m = realloc( *m , new_rows * sizeof (long long *));
		*row = new_rows;
	}
}

void matrix_free_long(long long ***m, int row, int col){
//Deallocazione di una matrice di tipo long long con dimensioni indicate.	
	for (int i=0; i<row; i++)      
		free((*m)[i]);
	free(*m);	
}

void matrix_alloc_int(int ***m, int row, int col){
//Allocazione di una matrice di tipo int con dimensioni indicate.	
	*m = malloc(row * sizeof (int *) );
	if( *m != NULL )
		for (int i=0; i<row; i++)
			(*m)[i] = calloc(col , sizeof (int) );	
}

void matrix_free_int(int ***m, int row, int col){
//Deallocazione di una matrice di tipo int con dimensioni indicate.	
	for (int i=0; i<row; i++)      
		free((*m)[i]);
	free(*m);	
}

//copia la matrice m1 nella matrice m2
void matrix_cpy(long long **m1, int row, int col, long long **m2){
	
	int i,j;
	for(i=0;i<row;i++){
		for(j=0;j<col;j++){
			m2[i][j] = m1[i][j];
		}
	}

}


void matrix_alloc_long(long long ***m, int row, int col){
//Allocazione di una matrice di tipo int con dimensioni indicate.	
	*m = malloc(row * sizeof (long long *) );
	if( *m != NULL )
		for (int i=0; i<row; i++){
			(*m)[i] = calloc(col , sizeof (long) );	
			if ((*m)[i]==NULL)
			{
				break;
			}
		}
}

void matrix_realloc_long(long long ***m, int new_row, int new_col){
	int i;
	*m = realloc( *m , (new_row) * sizeof (long long *));
	for(i=0; i<new_row; i++){
		(*m)[i] = realloc((*m)[i] , new_col * sizeof (long long) );
	}
}

//aggiunge la riga r alla matrice m, r deve avere linghezza uguale al numero delle colonne di m
void add_row_to_matrix(long long ***m, int *row, int col, long long *r){
	
	int i;
	*m = realloc( *m , (*row+1) * sizeof (long long *));
	(*m)[*row] = malloc(col * sizeof (long long) );
	for(i=0; i<col; i++){
		(*m)[*row][i] = r[i];
	}
	*row = *row + 1;	

}

void append_matrix(long long ***m1, int *row1, int col1, long long **m2, int row2, int col2){
	int i=0;
	if( col1 == col2 ){ //se le matrici hanno lo stesso numero di colonne
		for(i=0; i<row2; i++){
			add_row_to_matrix(m1,row1,col1,m2[i]);
		}
	}
}



void append_and_free_matrix(long long ***m1, int *row1, int col1, long long ***m2, int *row2, int col2){
	int i=0;
	int tot = *row1 + *row2;
	if( col1 == col2 ){ //se le matrici hanno lo stesso numero di colonne
		*m1 = realloc(*m1, tot * sizeof(long long *));
		for (int r = *row1; r < tot; r++) {
			(*m1)[r] = malloc(col1 * sizeof (long long));
			for (int c = 0; c < col1; c++)
				(*m1)[r][c] = (*m2)[r-*row1][c];
			free((*m2)[r-*row1]);
		}
		*row1 = tot;
		*row2 = 0;
		free(*m2);
	}
}

//funzione di confronto tra la righa rowA con rowB, scorrendo le colonne da destra a sinistra
//restituisce 1 se rowA > rowB, -1 se rowB > rowA, 0 altrimenti. Compatibile con qsort_r
//meglio per l'uguaglianza
int compare_rows(const void *rowA, const void *rowB, void *columns) {

	long long *row1, *row2;
	int col;
	
	col = *((int *) columns);
	row1 = *((long long **) rowA);
	row2 = *((long long **) rowB);
	
	for (int i = col-1; i >= 0; i--) {
		if (row1[i] > row2[i])
			return 1;
		else if (row1[i] < row2[i])
			return -1;
	}
	
	return 0;
}


//funzione di confronto tra la righa rowA con rowB, sommando il contenuto di ogni colonna
//restituisce 1 se rowA > rowB, -1 se rowB > rowA, 0 altrimenti. Compatibile con qsort_r
//meglio per il sorting
int compare_rows2(const void *rowA, const void *rowB, void *columns) {

	long long *row1, *row2;
	int col;
	int s1 = 0, s2 = 0;
	
	col = *((int *) columns);
	row1 = *((long long **) rowA);
	row2 = *((long long **) rowB);
	
	for (int i = col-1; i >= 0; i--) {
		s1 += (row1[i]*i);
		s2 += (row2[i]*i);
	}
	
	return s1-s2;
}


//tolgo da m1 le righe iniziali uguali a m2
void eliminate_equal_starting_rows(long long ***m1, int *row1, long long **m2, int row2, int col) {
		
		int eliminated_rows = 0, new_rows;
		long long **temp;
		
		//per ogni riga consecuitiva in m2 uguale a quella in m1
		while (eliminated_rows < row2 && !compare_rows(&m2[eliminated_rows], &(*m1)[eliminated_rows], &col))
			eliminated_rows++;
		
		//nuovo numero di righe
		new_rows = *row1 - eliminated_rows;	
		//alloco la nuova matrice e copio le righe non eliminate
		matrix_alloc_long(&temp, new_rows, col);
		for (int i = 0; i < new_rows; i++) {
			for (int j = 0; j < col; j++) {
				temp[i][j] = (*m1)[i+eliminated_rows][j];
			}
		}
		
		//elimino la vecchia matrice
		matrix_free_long(m1, *row1, col);
		
		*m1 = temp;
		*row1 = new_rows;
		
		return;
}

//elimina da m la righa di indice index
void delete_row(long long ***m, int *row, int col, int index) {
	
	//l'indice deve essere minore al numero di righe
	if (index >= *row)
		return;
	
	long long **temp;
	matrix_alloc_long(&temp, (*row)-1, col);
	
	//copio la matrice saltando la riga eliminata
	for (int r = 0; r < *row-1; r++)
		for (int c = 0; c < col; c++)
			if (r < index)
				temp[r][c] = (*m)[r][c];
			else
				temp[r][c] = (*m)[r+1][c];
			
	matrix_free_long(m, *row, col);
	
	*m = temp;
	*row = *row - 1;
	return;
}

//elimina da m1 le righe uguali a quelle di m2
void eliminate_equal_rows(long long ***m1, int *row1, long long **m2, int row2, int col) {
	
	//confronto ogni riga di m2 con ogni riga di m1
	for (int r = 0; r < row2; r++)
		for (int i = 0; i < *row1; i++)
			//se sono uguali elimino la riga da m1
			if (!compare_rows(&(*m1)[i], &m2[r], &col)) {
				delete_row(m1, row1, col, i);
				//avendone tolta una le righe sono diminuite, devo diminuire anche l'indice
				i--;
				//per risparmiare tempo assumo che tutte le righe di m1 siano diverse tra loro
				//possiamo farlo perchè usiamo questa funzione dopo gauss
				break;
			}
	return;
}

