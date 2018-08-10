#ifndef SCAN_H_  /* Include guard */
#define SCAN_H_

void allocation(int **matrix, int *row, int *col, int *numero_variabili, char **variabili, int *tipo_ordinamento, int *modulo, int *max_degree, FILE *input_file);

int parse(int num_var, char *vet, int **m, int row, int **vet_grd, int len, int module, int (*ord) (const void *, const void *, void*), FILE *input_file);

int parse_mon(char * mon, int len,int * val, int num_var, char *vet, int *grade, int module);

//associa ad *ord la funzione di ordinamento scelta
int order(int (**ord) (const void *, const void *, void*), int n);

#endif //SCAN_H_

