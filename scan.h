#ifndef SCAN_H_  /* Include guard */
#define SCAN_H_


void allocation(long long ***m, int *row, int *col, int *num_var, char **v, int *n, long long *module, int *max_degree, FILE *input_file);

int parse(int num_var, char *vet, long long **m, int row, int **vet_grd, int len, long long module, int (*ord) (const void *, const void *, void*), FILE *input_file);

int parse_mon(char * mon, int len,long long * val, int num_var, char *vet, int *grade, long long module);

//associa ad *ord la funzione di ordinamento scelta
int order(int (**ord) (const void *, const void *, void*), int n);

#endif //SCAN_H_

