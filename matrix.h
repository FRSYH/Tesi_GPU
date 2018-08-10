#ifndef MATRIX_H_  /* Include guard */
#define MATRIX_H_

#include <stdio.h>

void swap_rows(long long **m, int row, int col, int j, int i);  // scambia tra di loro due righe della matrice

void print_matrix (long long **m, int row, int col, FILE *output_file); // stampa la matrice

//copia il vettore vet2 in vet1, entrambi di lunghezza len
void vctcpy(int *vet1, int const *vet2, int len);

void matrix_free_long(long long ***m, int row, int col);

void matrix_alloc_int(int ***m, int row, int col);

void matrix_free_int(int ***m, int row, int col);

void matrix_cpy(long long **m1, int row, int col, long long **m2);

void matrix_alloc_long(long long ***m, int row, int col);

void matrix_realloc_long(long long ***m, int new_row, int new_col);

void add_row_to_matrix(long long ***m, int *row, int col, long long *r);

//compute the number of null rows (rows full of 0)
int null_rows(long long **m, int row, int col);

//eliminate the matrix null rows (reallocation - resize)
void eliminate_null_rows(long long ***m, int *row, int col);

void append_matrix(long long ***m1, int *row1, int col1, long long **m2, int row2, int col2);

void append_and_free_matrix(long long ***m1, int *row1, int col1, long long ***m2, int *row2, int col2);

//funzione di confronto tra la righa rowA con rowB, scorrendo le colonne da destra a sinistra
//restituisce 1 se rowA > rowB, -1 se rowB > rowA, 0 altrimenti. Compatibile con qsort_r
int compare_rows(const void *rowA, const void *rowB, void *columns);

int compare_rows2(const void *rowA, const void *rowB, void *columns);

//tolgo da m1 le righe iniziali uguali a m2
void eliminate_equal_starting_rows(long long ***m1, int *row1, long long **m2, int row2, int col);

//elimina da m la righa di indice index
void delete_row(long long ***m, int *row, int col, int index);

//elimina da m1 le righe uguali a quelle di m2
void eliminate_equal_rows(long long ***m1, int *row1, long long **m2, int row2, int col);


#endif //MATRIX_H_
