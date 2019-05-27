#ifndef LINALG_H_  /* Include guard */
#define LINALG_H_


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



#endif //LINALG_H_
