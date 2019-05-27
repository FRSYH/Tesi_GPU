#include <stdlib.h>

//n mod p 
//Riduzione di n in modulo p.
long long mod(long long n, long long p) {
	long long v = n, x = 0;

	if (v >= p) {
		v = n % p;
	}
	else {
		if (v < 0) {
			x = n / p;
			v = n - (x*p);
			v += p;
		}
	}
	return v;
}


//calcola il fattoriale di n (se n è negativo return -1)
long long factorial(int n) {
	long long k;

	if (n<0) //se n è negativo non esiste il fattoriale
	{
		return -1; //codice di errore
	}
	else { //altrimenti calcolo il fattoriale

		if (n == 0 || n == 1) {
			return 1;
		}
		else {
			k = 1;
			for (int i = 2; i <= n; i++) {
				k *= i;
			}
			return k;
		}
	}
}

//restituisce il numero di possibili monomi con n variabili e grado = m
int combination(int n, int m) {

	long long num, den;
	//calcolo {(m+n-1)! / m!*(n-1)!}

	//se n>=m semplificato a {(j+n-1)*(j+n-2)* ... *(n) / j!}
	if (n >= m) {
		num = 1;
		for (int k = m; k > 0; k--)
			num = num * (n + k - 1);
		den = factorial(m);
	}
	//se m>n semplificato a {(j+n-1)*(j+n-2)* ... *(j) / (n-1)!}
	else {
		num = 1;
		for (int k = n; k > 1; k--)
			num = num * (m + k - 1);
		den = factorial(n - 1);
	}
	return (int)(num / den);
}

//restituisce il numero di tutti i possibili monomi con n variabili e grado <= m
int monomial_combinations(int n, int m) {

	int result = 0;
	//result = Sommatoria (per j da 1 a m) {(j+n-1)! / j!*(n-1)!}
	for (int j = 0; j <= m; j++)
		result += combination(n, j);
	return  result;
}


//confronta due monomi di *arg variabili secondo l'ordinamento grevlex
//restituisce un intero positivo se monom1 > monom2, zero se sono uguali, uno negativo altrimenti
//i monomi sono sempre rappresentati come array di lunghezza pari al numero delle variabili
//sono fatti diversi cast perchè il tipo degli argomenti sia compatibile con qsort_r
int grevlex_comparison(const void *monom1, const void *monom2, void *arg) {

	int degree1 = 0, degree2 = 0, n, *mon1, *mon2;
	n = *((int *)arg);
	mon1 = *((int **)monom1);
	mon2 = *((int **)monom2);

	//calcolo i gradi dei monomi
	for (int v = 0; v < n; v++) {
		degree1 += mon1[v];
		degree2 += mon2[v];
	}
	if (degree1 > degree2)
		return 1;
	else if (degree1 < degree2)
		return -1;
	//se il grado è uguale guardo l'utlima cifra non nulla
	//del array risultante dalla sottrazione dei monomi
	else {
		int *temp = (int *)malloc(n * sizeof(int));
		int result;
		for (int v = 0; v < n; v++)
			temp[v] = mon1[v] - mon2[v];
		for (int v = (n - 1); v >= 0; v--) {
			if (temp[v] != 0) {
				result = -temp[v];
				free(temp);
				//per evitare di fare free due volte sul  puntatore lo setto a NULL dopo la free
				temp = NULL;
				return result;
			}
		}
		free(temp);
	}
	return 0;
}

//confronta due monomi di *arg variabili secondo l'ordinamento grevlex
//restituisce un intero positivo se monom1 > monom2, zero se sono uguali, uno negativo altrimenti
//i monomi sono sempre rappresentati come array di lunghezza pari al numero delle variabili
//sono fatti diversi cast perchè il tipo degli argomenti sia compatibile con qsort_r
int grevlex_comparison_mcvs(void *arg, const void *monom1, const void *monom2) {

	int degree1 = 0, degree2 = 0, n, *mon1, *mon2;
	n = *((int *)arg);
	mon1 = *((int **)monom1);
	mon2 = *((int **)monom2);

	//calcolo i gradi dei monomi
	for (int v = 0; v < n; v++) {
		degree1 += mon1[v];
		degree2 += mon2[v];
	}
	if (degree1 > degree2)
		return 1;
	else if (degree1 < degree2)
		return -1;
	//se il grado è uguale guardo l'utlima cifra non nulla
	//del array risultante dalla sottrazione dei monomi
	else {
		int *temp = (int *)malloc(n * sizeof(int));
		int result;
		for (int v = 0; v < n; v++)
			temp[v] = mon1[v] - mon2[v];
		for (int v = (n - 1); v >= 0; v--) {
			if (temp[v] != 0) {
				result = -temp[v];
				free(temp);
				//per evitare di fare free due volte sul  puntatore lo setto a NULL dopo la free
				temp = NULL;
				return result;
			}
		}
		free(temp);
	}
	return 0;
}

//https://git.devuan.org/jaretcantu/eudev/commit/a9e12476ed32256690eb801099c41526834b6390
//mancante nella stdlib, controparte di qsort_r
//effettua una ricerca binaria di key nell'array base di lunghezza nmemb i cui elementi
//hanno dimensione size, e restituisce un puntatore all'elemento uguale a key se c'è, altrimenti NULL.
//compar è la funzione di ordinamento con cui viene confrontato key con base
//arg è il terzo argomento di compar
void *bsearch_r(const void *key, const void *base, size_t nmemb, size_t size,
	int(*compar) (void *, const void *, const void *),
	void *arg) {
	size_t l, u, idx;
	const void *p;
	int comparison;

	l = 0;
	u = nmemb;
	while (l < u) {
		idx = (l + u) / 2;
		p = (void *)(((const char *)base) + (idx * size));
		comparison = compar(arg, key, p);
		if (comparison < 0)
			u = idx;
		else if (comparison > 0)
			l = idx + 1;
		else
			return (void *)p;
	}
	return NULL;
}