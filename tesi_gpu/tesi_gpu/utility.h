#pragma once
#ifndef UTILITY_H_  /* Include guard */
#define UTILITY_H_

struct map_row {
	int len;
	int *col;
};

struct map {
	int len;
	struct map_row *row;
};

void monomial_computation_rec(int n, int m, int **vet, int turn, int *monomial, int *pos);

int **monomial_computation(int n, int m, int len);

void setup_struct_map(struct map *map, int **monomi, int len, int n, int m, int(*compar) (void*, const void *, const void *));

int grado_monomio(int posizione, int **vet, int num_var);

void matrix_degree(int *m, int row, int col, int *m_deg, int **vet, int num_var);

int target_degree(int *v, int max_degree);

void init_degree_vector(int *degree, int num_var, int max_degree);

#endif //UTILITY_H_