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

#endif //UTILITY_H_