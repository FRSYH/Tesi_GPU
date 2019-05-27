#ifndef LINALG_H_  /* Include guard */
#define LINALG_H_
//n mod p
long long mod(long long n, long long p);

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

int grevlex_comparison_mcvs(void *arg, const void *monom1, const void *monom2);

//calcola il fattoriale di n
long long factorial(int n);

//mancante nella stdlib, controparte di qsort_r
void *bsearch_r(const void *key, const void *base, size_t nmemb, size_t size,
	int(*compar) (void *, const void *, const void *),
	void *arg);

#endif //LINALG_H_
