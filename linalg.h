#ifndef LINALG_H_  /* Include guard */
#define LINALG_H_

//#include <gmp.h>

//n mod p
long long mod(long long n, long long p);

// a + b mod p
long long add_mod(long long a, long long b, long long p);

// a - b mod p
long long sub_mod(long long a, long long b, long long p);

// a * b mod p
long long mul_mod(long long a, long long b, long long p);

// inverso moltiplicativo in mod p
long long invers(long long n, long long p);

// riduzione di gauss della matrice m
void magma_gauss(long long **m, int row, int col, int modulo);

// riduzione di gauss della matrice m
void gauss(long long **m, int row, int col, int modulo, int start, int *v);


void riduzione(long long **m, int row, int col, int riga_pivot, int j, int module);

//restituisce il numero di possibili monomi con n variabili e grad = m
int combination(int n, int m);

/*
//restituisce il numero di possibili monomi con n variabili e grad = m
int gmp_combination(int n, int m);
*/

//restituisce il numero di tutti i possibili monomi con n variabili e grado <=m
int monomial_combinations(int n, int m); 

//confronta due monomi di *arg variabili secondo l'ordinamento grevlex
//restituisce un intero positivo se monom1 > monom2, zero se sono uguali, uno negativo altrimenti
int grevlex_comparison(const void *mon1, const void *mon2, void *arg);

//calcola il fattoriale di n
long long factorial(int n);

/*
//calcola il fattoriale di n
void gmp_factorial(mpz_t result, int n);
*/

//mancante nella stdlib, controparte di qsort_r
void *bsearch_r(const void *key, const void *base, size_t nmemb, size_t size,
                 int (*compar) (const void *, const void *, void *),
                 void *arg);

#endif //LINALG_H_
