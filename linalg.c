#include <stdlib.h>
#include "matrix.h"
#include "linalg.h"
#include <stdio.h>
/*
#include <math.h>
#include <gmp.h>
*/

//n mod p 
//Riduzione di n in modulo p.
long long mod(long long n, long long p){
	long long v = n,x =0;

	if( v >= p ){
		v = n%p;
	}else{
		if( v < 0 ){
			x = n/p;
			v = n-(x*p);
			v += p;
		}
	}
	return v;
}

//inverso moltiplicativo di n in modulo p (con p primo).
long long invers(long long n, long long p){
	long long b0 = p, t, q;
	long long x0 = 0, x1 = 1;
	if (p == 1) return 1;
	while (n > 1) {
		q = n / p;
		t = p, p = (n % p), n = t;
		t = x0, x0 = x1 - q * x0, x1 = t;
	}
	if (x1 < 0) x1 += b0;
	return x1;	
}


// a + b mod p
//sommatoria di a e b in modulo p
long long add_mod(long long a, long long b, long long p){
	return mod((a+b),p);
}

// a - b mod p
//sottrazione di a e b in modulo p
long long sub_mod(long long a, long long b, long long p){
	return mod((a-b),p);
}

// a * b mod p
//prodotto di a e b in modulo p
long long mul_mod(long long a, long long b, long long p){
	return mod((a*b),p);
}


//restituisce il numero di tutti i possibili monomi con n variabili e grado <= m
int monomial_combinations(int n, int m) {

	int result = 0;
	//result = Sommatoria (per j da 1 a m) {(j+n-1)! / j!*(n-1)!}
	for (int j = 0; j <= m; j++)
		result += combination(n, j);
	return  result;
}

//Calcola la riduzione di Gauss della matrice m (matrice di grandezza row e col).
//La riduzione è calcolata sulla triangolare superiore sinistra.
void magma_gauss(long long **m, int row, int col, int modulo){
	
	int pivot_riga, pivot_colonna,righe_trovate,j;

	righe_trovate = -1;
	for( j=0; j<col; j++){
		pivot_colonna = col-(j+1);
		for( pivot_riga=(righe_trovate+1); pivot_riga<row; pivot_riga++ ){
			if( m[pivot_riga][pivot_colonna] != 0 ){
				riduzione(m,row,col,pivot_riga,pivot_colonna, modulo);
				righe_trovate++;
				swap_rows(m,row,col,righe_trovate,pivot_riga);
				break;
			}
		}
	}
}


/*  effettua la riduzione di gauss della matrice m in place
    parametro start viene utilizzato per ottimizzare l'agloritmo:
	- 0 effettua la normale riduzione di gauss
	- n (con n > 0) effettua una riduzione ottimizzata saltando le prime n righe nella fase di riduzione delle righe sottostanti.

	Questa ottimizzazione è possibile solo se le prime n righe sono risultato di una riduzione di gauss precedente, quindi risulta inutile ripetere l'operazione di riduzione su queste righe.
	Il salto delle prime n righe avviene solo se nella corrente riduzione non si è effettuato swap di righe. Quando le righe trovate superano start si ritorna alla normale riduzione.
	

*/
void gauss(long long **m, int row, int col, int module, int start, int *v){
	
	int pivot_riga = 0,r = 0,righe_trovate = 0,i,k;
	long long s,inv,a;
	int st,flag=0,invarianti=0,flag2=0,tmp;

	if( start == 0 ){
		flag = 1;
	}else{
		st = start;
	}

	for(int pivot_colonna = col-1; pivot_colonna >= 0; pivot_colonna-- ){
		r = righe_trovate;
		while( r < row && m[r][pivot_colonna] == 0 ){
			r++;
			
		}
		// ho trovato la prima riga con elemento non nullo in posizione r e pivot_colonna oppure non esiste nessuna riga con elemento non nullo in posizione pivot_colonna
		
		if( r < row ){ //significa che ho trovato un valore non nullo
			if( r != righe_trovate ){
				swap_rows(m,row,col,righe_trovate,r); //sposto la riga appena trovata nella posizone corretta
				flag = 1;
				if( v != NULL ){
					tmp = v[righe_trovate];
					v[righe_trovate] = v[r];
					v[r] = tmp;
				}				
			}			
			pivot_riga = righe_trovate;
			righe_trovate++;

			if( flag == 1 ){  			//riprendo la normale riduzione di gauss
				st = righe_trovate;
			}else{

				if( st < righe_trovate ){  //se sono nella modalitá ottimizzazione e supero le prime n righe allora ritorno alla normale riduzione
					flag = 1;
					st = righe_trovate;
				}
			}

			#pragma omp parallel for private(i,inv,s,k,a)	
			for( i = st; i < row; i++ ){
				if( m[i][pivot_colonna] != 0 ){
					inv = invers(m[pivot_riga][pivot_colonna],module);		//inverso dell´ elemento in m[r][pivot_colonna]
					s = mul_mod(inv,m[i][pivot_colonna],module);						
					//#pragma omp parallel for private (k,a) shared (m) schedule(dynamic)
					for( k = 0; k < pivot_colonna+1; k++ ){
						a = mul_mod(s,m[pivot_riga][k],module);
						m[i][k] = sub_mod(m[i][k],a,module);

					}
				}
			}
		}
	}
}



//Calcola la riduzione di Gauss di una singola riga della matrice m.
void riduzione(long long **m, int row, int col, int riga_pivot, int j, int module){
	
	int r,k;
	long long s,inv,a;
	#pragma omp parallel for private(inv,s,k,a)
	for( r=riga_pivot+1; r<row; r++ ){

		if( m[r][j] != 0 ){
			inv = invers(m[riga_pivot][j],module);			//calcola l'inverso moltiplicativo di m[riga_pivot][j] nel campo indicato
			s = mul_mod(inv,m[r][j],module);
			//#pragma omp parallel for private(a)
			for( k=0; k<col; k++ ){
				a = mul_mod(s,m[riga_pivot][k],module);
				m[r][k] = sub_mod(m[r][k],a,module);

			}
			
		}
	}
}


//confronta due monomi di *arg variabili secondo l'ordinamento grevlex
//restituisce un intero positivo se monom1 > monom2, zero se sono uguali, uno negativo altrimenti
//i monomi sono sempre rappresentati come array di lunghezza pari al numero delle variabili
//sono fatti diversi cast perchè il tipo degli argomenti sia compatibile con qsort_r
int grevlex_comparison(const void *monom1, const void *monom2, void *arg) {

	int degree1 = 0, degree2 = 0, n, *mon1, *mon2;
	n = *((int *) arg);
	mon1 = *((int **) monom1);
	mon2 = *((int **) monom2);
	
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
		int *temp = malloc(n * sizeof(int));
		int result;		
		for (int v = 0; v < n; v++)
			temp[v] = mon1[v] - mon2[v];
		for (int v = (n-1); v >= 0; v--) {
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


//calcola il fattoriale di n (se n è negativo return -1)
long long factorial(int n){
	long long k;

	if (n<0) //se n è negativo non esiste il fattoriale
	{
		return -1; //codice di errore
	}else{ //altrimenti calcolo il fattoriale

		if( n==0 || n==1 ){
			return 1;
		}else{
			k=1;
			for (int i = 2; i <= n; i++){
				k *= i;	
			}
			return k;
		}
	}
}

/*
void gmp_factorial(mpz_t result, int n){

	mpz_t x;
	mpz_init(x);
	if (n<0) //se n è negativo non esiste il fattoriale
	{
		mpz_set_si(result,-1);

	}else{ //altrimenti calcolo il fattoriale

		if( n==0 || n==1 ){
			mpz_set_si(result,1);
		}else{
			
			mpz_set_si(result,1);
			for (int i = 2; i <= n; i++){
				mpz_set_si(x,i);
				mpz_mul(result,result,x);
			}
		}
	}
}
*/ 

//restituisce il numero di possibili monomi con n variabili e grado = m
int combination(int n, int m){

	long long num, den;
	//calcolo {(m+n-1)! / m!*(n-1)!}

	//se n>=m semplificato a {(j+n-1)*(j+n-2)* ... *(n) / j!}
	if (n >= m) {
		num = 1;
		for (int k = m; k > 0; k--)
			num = num * (n+k-1);
		den = factorial(m);
	}
	//se m>n semplificato a {(j+n-1)*(j+n-2)* ... *(j) / (n-1)!}
	else {
		num = 1;
		for (int k = n; k > 1; k--)
			num = num * (m+k-1);
		den = factorial(n-1);
	}
	return (num/den);		
}


/*
//restituisce il numero di possibili monomi con n variabili e grado = m
int gmp_combination(int n, int m){

	mpz_t num, den, x, result;
	mpz_init(num);
	mpz_init(den);
	mpz_init(x);
	mpz_init(result);
	
	//calcolo {(m+n-1)! / m!*(n-1)!}

	//se n>=m semplificato a {(j+n-1)*(j+n-2)* ... *(n) / j!}
	if (n >= m) {
		
		mpz_set_ui(num,1); 
		for (int k = m; k > 0; k--){
			mpz_set_si(x,n+k-1);	
			mpz_mul(num,num,x);
		}	
		
		gmp_factorial(den,m);
	}
	//se m>n semplificato a {(j+n-1)*(j+n-2)* ... *(j) / (n-1)!}
	else {
		mpz_set_ui(num,1);
		for (int k = n; k > 1; k--){
			mpz_set_si(x,m+k-1);	
			mpz_mul(num,num,x);			
		}
		
		gmp_factorial(den,n-1);
	}
	
	mpz_divexact(result,num,den);	// result = num/den
	if( mpz_fits_sint_p(result) != 0 ){
		int r = 0;
		mpz_export(&r, 0, -1, sizeof result, 0, 0, result);
		return r;
		//printf("Trovato risultato che entra in int: %d\n", r);

	}else{
		printf("Overflow gmp_combination\n");
		return 0;
	}
}
*/

//https://git.devuan.org/jaretcantu/eudev/commit/a9e12476ed32256690eb801099c41526834b6390
//mancante nella stdlib, controparte di qsort_r
//effettua una ricerca binaria di key nell'array base di lunghezza nmemb i cui elementi
//hanno dimensione size, e restituisce un puntatore all'elemento uguale a key se c'è, altrimenti NULL.
//compar è la funzione di ordinamento con cui viene confrontato key con base
//arg è il terzo argomento di compar
void *bsearch_r(const void *key, const void *base, size_t nmemb, size_t size,
                 int (*compar) (const void *, const void *, void *),
                 void *arg) {
	size_t l, u, idx;
	const void *p;
	int comparison;

	l = 0;
	u = nmemb;
	while (l < u) {
		idx = (l + u) / 2;
		p = (void *)(((const char *) base) + (idx * size));
		comparison = compar(key, p, arg);
		if (comparison < 0)
			u = idx;
		else if (comparison > 0)
			l = idx + 1;
		else
			return (void *)p;
	}
	return NULL;
}

