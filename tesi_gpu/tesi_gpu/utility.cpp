#include <stdlib.h>
#include "matrix.h"
#include "scan.h"
#include "linalg.h"
#include "utility.h"

/*funzione ricorsiva che calcola tutti i possibili monomi con n variabili e grado <= m
e li inserisce nell'array vet. I monomi sono rappresentati come array di interi dove
il valore di ogni posizione rappresenta il grado della variabile in quella posizione.
Esempio: n=3, x^2*y*z = [2,1,1].
L'array vet deve essere già allocato correttamente. Gli altri parametri sono necessari
per la struttura ricorsiva della funzione e alla prima chiamata devono essere:
- turn = 0, rappresenta la posizione della variabile nel monomio
- monomial = array di interi di lunghezza n già allocato e usato per calcolare i vari monomi
- *pos = 0 puntatore ad intero, rappresenta la prima posizione libera nell'array vet
*/

void monomial_computation_rec(int n, int m, int **vet, int turn, int *monomial, int *pos) {

	//per ogni variabile provo tutti i gradi da 0 a m
	for (int degree = 0; degree <= m; degree++) {
		//se questa è la prima variabile azzero il monomio
		if (turn == 0) {
			//azzero il monomio lasciando solo il grado della prima variabile
			monomial[0] = degree;
			for (int v = 1; v < n; v++)
				monomial[v] = 0;
		}
		//altrimenti le altre variabili aggiungo il proprio grado al monomio
		else
			monomial[turn] = degree;


		//ottengo il grado del monomio sommando i gradi delle variabili
		int sum = 0;
		for (int v = 0; v <= turn; v++)
			sum += monomial[v];
		//se il grado del monomio supera quello massimo non ha senso continuare a cercare
		//altri monomi partendo da questo, perchè tutti avranno grado maggiore o uguale
		if (sum > m)
			break;

		//se questa è l'ultima variabile copia il monomio nell'array vet and incrementa l'indice pos
		if (turn == (n - 1)) {
			vctcpy(vet[(*pos)], monomial, n);
			(*pos)++;
		}
		//altrimenti richiama se stessa cambiando la variabile (turn)
		else
			monomial_computation_rec(n, m, vet, turn + 1, monomial, pos);
	}

	return;
}

/*restituisce un array contenente tutti i len monomi con n variabili e grado <= m
len è il numero di possibili monomi con n variabili e grado <= m
i monomi sono array di interi di lunghezza n dove il valore di ogni posizione rappresenta
il grado della variabile in quella posizione. Esempio: n=3, x^2*y*z = [2,1,1]
len viene passato come argomento per evitare di ricalcolarlo internamente
*/
int **monomial_computation(int n, int m, int len) {

	int **vet, *monomial;

	//alloco la memoria per l'array
	matrix_alloc_int(&vet, len, n);

	//strutture di supporto necessarie per il calcolo
	monomial = (int *)malloc(n * sizeof(int));
	int pos = 0;

	//il calcolo è fatto dalla funzione ricorsiva correttemente parametrizzata
	monomial_computation_rec(n, m, vet, 0, monomial, &pos);

	free(monomial);

	return vet;
}

void setup_struct_map(struct map *map, int **monomi, int len, int n, int m, int(*compar) (void*, const void *, const void *)) {

	int sum, index = len;

	//	inizializzo la struttura map, la mappa ha len righe.
	map->len = len;
	map->row = (map_row *)malloc(map->len * sizeof(struct map_row));

	//per ogni monomio in vet
	int row, col, i, v;
	for (row = 0; row < len; row++) {
		index = 0;
		//dichiarati dentro per la parallelizzazione
		int *temp = (int *)malloc(n * sizeof(int));
		int *save = (int *)calloc(len, sizeof(int));
		//provo a moltiplicarlo con ogni monomio in vet
		for (col = 0; col < len; col++) {
			sum = 0;
			//eseguo il prodotto (sum è la somma dei gradi)
			for (v = 0; v < n; v++) {
				temp[v] = monomi[row][v] + monomi[col][v];
				sum += temp[v];
			}
			//se il grado del prodotto > grado massimo tutti i restanti prodotti
			//su quella riga sono > grado massimo
			if (sum > m) {

				//	a questo punto col è l'indice del primo elemento della mappa che non è possibile rappresentare, quindi la riga row ha solo col numero di celle e non len come prima.
				index = col;
				break;
			}
			//altrimenti cerco il prodotto in vet e metto l'indice in save
			else {
				save[col] = int((int **)(bsearch_r((void *)&temp, (void *)monomi, len, (sizeof(int*)), compar, &n)) - monomi);
			}
		}

		//	terminato il ciclo sulle colonne posso inizializzare la struttura perchè conosco tutti gli elementi da inserire	
		//  la riga attuale ha esattamente index elementi diversi da -1, quindi la riga avrà lunghezza pari a index precedentemente calcolato
		//  alloco la riga con un array da index elementi

		map->row[row].len = index;
		map->row[row].col = (int *)malloc(map->row[row].len * sizeof(int));
		//	a questo map devo copiare gli elementi generati dento alla struttura

		for (i = 0; i<map->row[row].len; i++)
			map->row[row].col[i] = save[i];

		free(temp);
		free(save);
	}
}

int grado_monomio(int posizione, int **vet, int num_var) {
	//Calcola il grado del monomio a partire dalla posizione occupata nel vettore (ordinato) delle posizioni rispetto l'ordinamento scelto.
	//(la posizione occupata deve essere corretta).
	int i, grado;
	grado = 0;
	for (i = 0; i<num_var; i++) {
		grado += vet[posizione][i];
	}
	return grado;
}

void matrix_degree(int *m, int row, int col, int *m_deg, int **vet, int num_var) {
	//m_deg è un vettore che ha lunghezza pari al grado massimo.
	//la funzione calcola i gradi dei polinomi presenti nella matrice.
	//Ogni cella del vettore m_deg rappresenta un grado, se esso compare nella matrice allora viene impostato a 1 o altrimenti.

	int i, j, last, grado, linear_index = 0;
	for (i = 0; i<row; i++) {
		for (j = col - 1; j>0; j--) {
			linear_index = i * col + j;
			if (m[linear_index] != 0) {
				last = j;           //posizione dell'ultimo coefficiente della riga
				break;
			}
		}
		grado = grado_monomio(last, vet, num_var);
		m_deg[grado] = 1;
	}
}

int target_degree(int *v, int max_degree) {
	//Controlla se il vettore v rappresenta la condizione di terminazione con gradi completi {1,2,3,...,max_degree}
	//Se la condizione è soddisfatta return 0 altrimenti -1

	int i, flag;
	flag = 0;
	for (i = 1; i<max_degree + 1; i++) {
		if (v[i] != 1) {
			flag = -1;
			break;
		}
	}
	return flag;
}

void init_degree_vector(int *degree, int num_var, int max_degree) {
	//inizializza il vettore degree con il numero di monomi di grado i-esimo <= del grado massimo
	int i, c;
	for (i = 0; i<max_degree + 1; i++) {
		c = combination(num_var, i);
		degree[i] = c;
	}
}