#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <time.h>
#include <omp.h>
#include <stdbool.h>
#include <time.h>
#include <cuda_runtime.h>

// compilazione nvcc gm.cu -o gm -w -Xcompiler " -openmp"
// nvcc gm.cu -o gm -w -Xcompiler " -openmp" -gencode arch=compute_61,code=sm_61 -lcudadevrt -rdc=true



//dichiarazione variabili globali
int max_degree = 0;
int module = 0;


struct map_row {
	int len;
	int *col;
};

struct map {
	int len;
	struct map_row *row;
};


//----------------------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------------------


void matrix_alloc_int(int ***m, int row, int col) {
	//Allocazione di una matrice di tipo int con dimensioni indicate.	
	*m = (int **)malloc(row * sizeof(int *));
	if (*m != NULL)
		for (int i = 0; i<row; i++)
			(*m)[i] = (int *)calloc(col, sizeof(int));
}

void matrix_free_int(int ***m, int row, int col) {
	//Deallocazione di una matrice di tipo int con dimensioni indicate.	
	for (int i = 0; i<row; i++)
		free((*m)[i]);
	free(*m);
}

//copia il vettore vet2 in vet1, entrambi di lunghezza len
void vctcpy(int *vet1, const int *vet2, int len) {
	for (int i = 0; i < len; i++)
		vet1[i] = vet2[i];
	return;
}

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
	matrix_alloc_int(&vet,len,n);

	//strutture di supporto necessarie per il calcolo
	monomial = (int *)malloc(n * sizeof(int));
	int pos = 0;

	//il calcolo è fatto dalla funzione ricorsiva correttemente parametrizzata
	monomial_computation_rec(n, m, vet, 0, monomial, &pos);

	free(monomial);

	return vet;
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
	return (num / den);
}

//restituisce il numero di tutti i possibili monomi con n variabili e grado <= m
int monomial_combinations(int n, int m) {

	int result = 0;
	//result = Sommatoria (per j da 1 a m) {(j+n-1)! / j!*(n-1)!}
	for (int j = 0; j <= m; j++)
		result += combination(n, j);
	return  result;
}

void allocation(int **matrix, int *row, int *col, int *numero_variabili, char **variabili, int *tipo_ordinamento, int *modulo, int *max_degree, FILE *input_file) {
	/*
	Legge da input le seguenti informazioni:
	- modulo dei coefficienti
	- grado massimo
	- numero dei polinomi di partenza
	- tipo di ordinamento
	- variabili utilizzate nei polinomi
	con queste informazioni alloca la matrice principale (matrice che conterrà i polinomi) e stabilisce il numero di variabili utilizzate.
	*/
	fscanf(input_file, "%d", modulo); //leggo il modulo
	fgetc(input_file);
	fscanf(input_file, "%d", max_degree); //leggo il grado massimo
	fgetc(input_file);
	fscanf(input_file, "%d", row);  //leggo numero dei polinomi di partenza
	fgetc(input_file);
	fscanf(input_file, "%d", tipo_ordinamento);  //leggo tipo di ordinamento
	fgetc(input_file);

	int i, j, k, pos_pol, num_pol;
	char c;

	i = 0;
	pos_pol = 0;
	*variabili = (char *)malloc(sizeof(char));
	c = fgetc(input_file);
	while (c != '\n') {
		(*variabili)[i] = c;
		i++;
		(*numero_variabili)++;
		*variabili = (char *)realloc(*variabili, (i + 1) * sizeof(char));
		c = fgetc(input_file);
	}

	*col = monomial_combinations(*numero_variabili, *max_degree);
	*matrix = (int *)calloc((*row) * (*col), sizeof(int));

}


void print_matrix(int *matrix, int row, int col, FILE *output_file) {
	for (int x = 0; x < row; x++) {
		for (int y = 0; y < col; y++) {
			fprintf(output_file, "%d ", matrix[ (x*col) + y]);
		}
		fprintf(output_file, "\n\n\n");
	}
	fprintf(output_file, "\n");
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

int order(int(**ord) (void*, const void *, const void *), int n) {
	//inizializza il puntatore ord alla funzione di ordinamento adeguata. Il numero n indica quale funzione scegliere.

	switch (n) {

	case 0:
		*ord = grevlex_comparison_mcvs;
		return 0;
		break;

	default:
		return -1;
		break;

	}
}


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



/* mon: stringa che rappresenta un monomio (non c'è carattere terminazione stringa)
* len: numero di caratteri in mon
* val: variabile in cui salvare il coefficiente del monomio
* num_var: numero di variabili nel sistema
* vet: vettore di caratteri in cui ogni carattere è una variabile (letto precedentemente da input)
* grade: vettore in cui salvare i gradi delle variabili secondo l'ordine di vet
* module: campo su cui è rappresentato il sistema
*/

int parse_mon(char * mon, int len, int * val, int num_var, char *vet, int *grade, int module) {

	//parsing prima del coefficiente
	int index = 0;
	//se il primo carattere è una lettera (variabile) il coefficiente è 1
	if (isalpha(mon[index]))
		*val = 1;
	//altrimenti leggo il coefficiente
	else {
		//se non è nè lettera nè cifra il formato input è sbagliato
		if (!isdigit(mon[index]))
			return -1;

		char *coefficient = (char *)malloc(sizeof(char));
		while (index < len && isdigit(mon[index])) {
			coefficient = (char *)realloc(coefficient, (index + 1) * sizeof(char));
			coefficient[index] = mon[index];
			index++;
		}
		//aggiungo il carattere di temrinazione
		coefficient = (char *)realloc(coefficient, (index + 1) * sizeof(char));
		coefficient[index] = '\0';
		//traduco il coefficiente in valore numerico e calcolo il modulo
		*val = mod(atoll(coefficient), module);
		free(coefficient);
	}

	//assumo grado zero di ogni variabile, aggiornato in seguito
	for (int k = 0; k < num_var; k++)
		grade[k] = 0;

	//parsing delle incognite
	char variable;
	int variable_degree;
	int variable_index;
	int exponent_index;
	char *exponent;

	while (index < len) {
		variable_index = -1;
		variable_degree = 0;

		//salto il moltiplicatore
		if (mon[index] == '*' || mon[index] == ' ')
			index++;
		//leggo la variabile
		if (index < len && isalpha(mon[index])) {
			variable = mon[index];

			//cerco la posizione della variabile in vet
			for (int i = 0; i < num_var; i++)
				if (vet[i] == mon[index]) {
					variable_index = i;
					//se è presente ha almeno grado 1
					variable_degree = 1;
					break;
				}
			//se non trovo la variabile in vet segnalo errore
			if (variable_index == -1)
				return -1;
			index++;
		}

		//se c'è il carattere dell'elevato leggo l'esponente
		if (index < len && mon[index] == '^') {
			index++;
			exponent_index = 0;
			//se non è una cifra segnalo errore
			if (index > len || !isdigit(mon[index]))
				return -1;
			exponent = (char *)malloc(sizeof(char));
			while (index < len && isdigit(mon[index])) {
				exponent = (char *)realloc(exponent, (exponent_index + 1) * sizeof(char));
				exponent[exponent_index] = mon[index];
				exponent_index++;
				index++;
			}
			//metto il carattere di terminazoine stringa
			exponent = (char *)realloc(exponent, (exponent_index + 1) * sizeof(char));
			exponent[exponent_index] = '\0';
			//leggo l'esponente
			variable_degree = atoi(exponent);
			free(exponent);
		}
		//se c'è la variabile	
		if (variable_index != -1)
			//metto in grade il grado della variabile nella posizione corretta
			grade[variable_index] = variable_degree;
	}
	return 0;
}



int parse(int numero_variabili, char *variabili, int *matrix, int row, int **monomi, int len, int module, int(*ord) (void*, const void *, const void *), FILE *input_file) {
	/*
	Esegue la lettura (parse) dei polinomi di partenza nel seguente modo.
	Si legge un monomio alla volta.
	Il monomio viene scomposta dalla funzione parse_mon.
	Si inserisce il coefficiente del monomio nella matrice principale (matrice dei coefficienti) nella posizione corretta.
	La posizione corretta è indicata da vet_grd.
	Si leggono tutti i monomi di tutti i polinomi di partenza.
	In caso di errore di formato nell'input la funzione si interrompe restituendo segnale di errore -1.
	*/

	int pos_pol = 0, i, col;
	char c, *mon;
	int cof = 0;
	c = fgetc(input_file);
	int linear_index = 0;
	int *grade;

	//finchè non termino il file o non ho terminato il numero di polinomi dichiarati
	while (c != EOF && pos_pol < row) {
		mon = (char *)malloc(sizeof(char));
		grade = (int *)calloc(numero_variabili, sizeof(int));
		i = 0;
		while (c != '+' && c != EOF && c != '\n') {
			mon = (char *)realloc(mon, (i + 1) * sizeof(char));
			mon[i] = c;
			i++;
			c = fgetc(input_file);
		}
		//se non ho salvato niente in mon (i = 0) non faccio il parsing
		if (i != 0 && parse_mon(mon, i, &cof, numero_variabili, variabili, grade, module) == -1) {
			return -1;
		}
		//inserire monomio in posizione corretta
	
		col = int((int **)(bsearch_r((void *)&grade, (void *)monomi, len, (sizeof(int*)), ord, &numero_variabili)) - monomi);
		linear_index = (pos_pol * len) + col;
		matrix[linear_index] = cof;
		if (c == '\n') {
			pos_pol++;
		}
		free(mon);
		free(grade);
		c = fgetc(input_file);
		cof = 0;
	}
	return 0;
}

int init_matrix(int *matrix, int row, int col, int **vet_grd, char *variabili, int num_var, int(*ord) (void*, const void *, const void *), FILE *input_file) {
	//Inizializza la matrice principale (dei coefficienti) con i coefficienti dei polinomi forniti come input.
	return parse(num_var, variabili, matrix, row, vet_grd, col, module, ord, input_file);
}

void setup_struct_map(struct map *map, int **monomi, int len, int n, int m, int (*compar) (void*, const void *, const void *)  ){
	
	int sum, index=len;

	//	inizializzo la struttura map, la mappa ha len righe.
	map->len = len;
	map->row = (map_row *)malloc( map->len * sizeof(struct map_row) );

	//per ogni monomio in vet
	int row, col, i, v;
	for (row = 0; row < len; row++){
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
			else{
				save[col] = (int **)(bsearch_r((void *) &temp, (void *) monomi, len, (sizeof(int*)), compar, &n)) - monomi;
			}
		}

		//	terminato il ciclo sulle colonne posso inizializzare la struttura perchè conosco tutti gli elementi da inserire	
		//  la riga attuale ha esattamente index elementi diversi da -1, quindi la riga avrà lunghezza pari a index precedentemente calcolato
		//  alloco la riga con un array da index elementi

		map->row[row].len = index;
		map->row[row].col = (int *)malloc( map->row[row].len * sizeof(int) );
		//	a questo map devo copiare gli elementi generati dento alla struttura

		for(i=0; i<map->row[row].len; i++)
			map->row[row].col[i] = save[i];
		
		free(temp);
		free(save);
	}
}

void init_degree_vector(int *degree, int num_var){
//inizializza il vettore degree con il numero di monomi di grado i-esimo <= del grado massimo
	int i,j,c;
	for(i=0; i<max_degree+1; i++){
		c = combination(num_var,i);
		degree[i] = c;
	}
}

int grado_monomio(int posizione, int **vet, int num_var){
//Calcola il grado del monomio a partire dalla posizione occupata nel vettore (ordinato) delle posizioni rispetto l'ordinamento scelto.
//(la posizione occupata deve essere corretta).
	int i,grado;
	grado = 0;
	for(i=0; i<num_var; i++){
		grado += vet[posizione][i];
	}
	return grado;
}

void matrix_degree(int *m, int row, int col, int *m_deg, int **vet, int num_var){
//m_deg è un vettore che ha lunghezza pari al grado massimo.
//la funzione calcola i gradi dei polinomi presenti nella matrice.
//Ogni cella del vettore m_deg rappresenta un grado, se esso compare nella matrice allora viene impostato a 1 o altrimenti.

	int i,j,last,grado, linear_index = 0;
	for(i=0; i<row; i++){
		for(j=col-1; j>0; j--){
			linear_index = i*col + j;
			if( m[linear_index] != 0 ){
				last = j;           //posizione dell'ultimo coefficiente della riga
				break;
			}
		}
		grado = grado_monomio(last, vet, num_var);
		m_deg[grado] = 1;
	}
}

void moltiplica_matrice(int **m, int *row, int col, struct map map, int * degree, int **vet, int num_var, int expansion_degree){
	
	int riga;
	int grado_massimo_riga, grado_massimo_monomio,i,j,last,new_row = 0;
	last = -1;
	int linear_index = 0;
	long long total_dim = 0;
	int *last_index = (int*)calloc(*row, sizeof(int));
	int *numero_polinomi = (int*)calloc(*row, sizeof(int));
	int numero_nuovi_polinomi = 0;

	for(riga=0; riga<*row; riga++){
		for(i=col-1; i>0; i--){
			linear_index = riga * col + i;
			if( (*m)[linear_index] != 0 ){  //(*m)[riga][i] != 0
				last = i;
				break;
			}
		}
		//risalgo al grado del monomio appena trovato
		//scorro la lista delle posizioni di inizio dei monomi con lo stesso grado
		
		last_index[riga] = last;

		if( last != -1 ){

			grado_massimo_riga = grado_monomio(last,vet,num_var);

			//calcolo il grado massimo che deve avere il monomio per cui moltiplicare
			grado_massimo_monomio = max_degree - grado_massimo_riga;
			// a questo punto conosco per quanti monomi devo moltiplicare e quindi
			// conosco il numero di righe che devo aggiungere alla matrice
			if( expansion_degree != 0 ){
				if( grado_massimo_monomio > expansion_degree ){
					grado_massimo_monomio = expansion_degree;
				}
			}

			for(i=1; i<(grado_massimo_monomio+1); i++){
				new_row += degree[i];
				numero_nuovi_polinomi += degree[i];
			}
			numero_polinomi[riga] = numero_nuovi_polinomi;
			numero_nuovi_polinomi = 0;
		}
	}
	
//--------------------------------------------------------------
	
	//printf("nuove righe %d, totale righe %d", new_row, (*row+new_row) );

	//ridimensionamento della matrice
	total_dim = (*row * col) + (new_row * col);
	*m = (int *)realloc( *m , total_dim * sizeof(int) );
	//azzeramento delle nuove righe
	for(i=*row; i<*row+new_row; i++){
		for(j=0; j<col; j++){
			(*m)[i*col+j] = 0;
		}
	}

	int len = *row;
	for(riga=0; riga<len; riga++){
		if( last_index[riga] != -1 ){
			for(i=1; i<(numero_polinomi[riga]+1); i++){     								//scorre tutti i monomi per i quali posso moltiplicare
				for(j=0; j<(last_index[riga]+1); j++){     								//scorre fino all'ultimo elemento della riga
					//(*m)[*row][ map.row[i].col[j] ] = (*m)[riga][j];  				//shift nella posizione corretta indicata dalla mappa
					linear_index = *row * col + map.row[i].col[j];
					(*m)[linear_index] = (*m)[riga*col+j];				
				}
				*row = *row + 1;											//aumento del conteggio delle righe
			}			
		}	
	}
	
	free(last_index);
	free(numero_polinomi);
}

void moltiplica_riga_forn(int **m, int *row, int col, int riga, struct map map, int * degree, int **vet, int num_var, int stop_degree){

	int grado_massimo_riga, grado_massimo_monomio,i,j,last,new_row;
	last = -1;
	int linear_index = 0;
	long long total_dim = 0;
	//cerco la posizione dell'ultimo coefficiente non nullo del polinomio rappresentato nella riga.
	for(i=col-1; i>0; i--){
		linear_index = riga * col + i;
		if( (*m)[linear_index] != 0 ){  //(*m)[riga][i] != 0
			last = i;
			break;
		}
	}
	//risalgo al grado del monomio appena trovato
	//scorro la lista delle posizioni di inizio dei monomi con lo stesso grado
	if( last != -1 ){

		grado_massimo_riga = grado_monomio(last,vet,num_var);

		//calcolo il grado massimo che deve avere il monomio per cui moltiplicare
		grado_massimo_monomio = max_degree - grado_massimo_riga;
		// a questo punto conosco per quanti monomi devo moltiplicare e quindi
		// conosco il numero di righe che devo aggiungere alla matrice
		new_row = 0;

		for(i=1; i<(grado_massimo_monomio+1); i++){
			new_row += degree[i];
		}

		total_dim = (*row * col) + (new_row * col);
		*m = (int *)realloc( *m , total_dim * sizeof(int) );
		//azzeramento delle nuove righe
		for(i=*row; i<*row+new_row; i++){
			for(j=0; j<col; j++){
				(*m)[i*col+j] = 0;
			}
		}

		for(i=1; i<(new_row+1); i++){     								//scorre tutti i monomi per i quali posso moltiplicare
			for(j=0; j<(last+1); j++){     								//scorre fino all'ultimo elemento della riga
				//(*m)[*row][ map.row[i].col[j] ] = (*m)[riga][j];  				//shift nella posizione corretta indicata dalla mappa
				linear_index = *row * col + map.row[i].col[j];
				(*m)[linear_index] = (*m)[riga*col+j];				
			}
			*row = *row + 1;											//aumento del conteggio delle righe
		}
	}

}

//Scambio di due righe della matrice m.
void swap_rows(int *m, int row, int col, int j, int i){
	
	int k;
	long long tmp;
	if( j!=i ){
		for(k=0;k<col;k++){
			tmp = m[i*col+k];			//m[i][k];
			m[i*col+k] = m[j*col+k];	//m[i][k] = m[j][k];
			m[j*col+k] = tmp;			//m[j][k] = tmp;
		}
	}
}


//Scambio di due righe della matrice m.
__device__ void swap_rows_GPU(int *m, int row, int col, int j, int i){
	
	int k;
	long long tmp;
	if( j!=i ){
		for(k=0;k<col;k++){
			tmp = m[i*col+k];			//m[i][k];
			m[i*col+k] = m[j*col+k];	//m[i][k] = m[j][k];
			m[j*col+k] = tmp;			//m[j][k] = tmp;
		}
	}
}

int mod_long(long long n, int p) {
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
	int r = v;
	return r;
}


//n mod p 
//Riduzione di n in modulo p.
__device__ int mod_long_GPU(long long n, long long p) {
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
	int r = v;
	return r;
}


//n mod p 
//Riduzione di n in modulo p.
__device__ int mod_GPU(int n, int p) {
	int v = n, x = 0;

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

//inverso moltiplicativo di n in modulo p (con p primo).
int invers(int n, int p){
	int b0 = p, t, q;
	int x0 = 0, x1 = 1;
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
int add_mod(int a, int b, int p){
	return mod((a+b),p);
}

// a - b mod p
//sottrazione di a e b in modulo p
int sub_mod(int a, int b, int p){
	long long aa,bb;
	aa = a;
	bb = b;
	return mod_long((aa-bb),p);
}

// a * b mod p
//prodotto di a e b in modulo p
int mul_mod(int a, int b, int p){
	long long aa,bb;
	aa = a;
	bb = b;
	return mod_long((aa*bb),p);
}


//inverso moltiplicativo di n in modulo p (con p primo).
__device__ int invers_GPU(int n, int p){
	int b0 = p, t, q;
	int x0 = 0, x1 = 1;
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
__device__ int add_mod_GPU(int a, int b, int p){
	return mod_GPU((a+b),p);
}

// a - b mod p
//sottrazione di a e b in modulo p
__device__ int sub_mod_GPU(int a, int b, int p){
	long long aa,bb;
	aa = a;
	bb = b;
	return mod_long_GPU((aa-bb),p);
}

// a * b mod p
//prodotto di a e b in modulo p
__device__ int mul_mod_GPU(int a, int b, int p){
	long long aa,bb;
	aa = a;
	bb = b;
	return mod_long_GPU((aa*bb),p);
}



#define gpuErrchk(ans) { gpuAssert((ans), __FILE__, __LINE__); }
inline void gpuAssert(cudaError_t code, const char *file, int line, bool abort=true)
{
   if (code != cudaSuccess) 
   {
      fprintf(stderr,"GPUassert: %s %s %d\n", cudaGetErrorString(code), file, line);
      if (abort) exit(code);
   }
}


/*
	Ottimizzazioni
		- branch divergence

	Correttezza
		- kernel celle (calcoli errati)

	Aggiunte
		- risoluzione e stampa soluzioni delle incognite

	Test
		- performance al variare della dimensione di grid e block

*/


__global__ void kernel_riduzione_riga(int *matrix, int row, int col, int dim, int module, int start, int pivot_colonna, int inv, int pivot_riga){

	int a = 0, s = 0;
	int last_row = row - 1;
	//int row_index =  last_row - (blockDim.x * blockIdx.x + threadIdx.x); 
	int row_index = (pivot_riga + 1) + (blockDim.x * blockIdx.x + threadIdx.x);
	if(row_index >= start && row_index < row){
		
		int row_linear_index = row_index * col + pivot_colonna;
		if( matrix[row_linear_index] != 0 ){					
			s = mul_mod_GPU(inv,matrix[row_linear_index],module);						
			for(int k = 0; k < pivot_colonna+1; k++ ){
				a = mul_mod_GPU(s,matrix[pivot_riga*col+k],module);
				matrix[row_index*col+k] = sub_mod_GPU(matrix[row_index*col+k],a,module);		
			}

		}
	}
}


__global__ void kernel_riduzione_cella(int *matrix, int row, int col, int dim, int module, int inv, int pivot_colonna, int pivot_riga, int start){

	int last_row = row - 1;
	int starting_row = pivot_riga + 1;
	int idx = blockIdx.x * blockDim.x + threadIdx.x;
	int idy = blockIdx.y * blockDim.y + threadIdx.y + starting_row; 

	if( idx < pivot_colonna  && idy < row && idy > pivot_riga){		//fermo i thread prima di pivot_colonna per impedire di sovrascrivere il dato necessario per s
		int div = matrix[idy*col+pivot_colonna];
		if( div != 0 ){
			int s = mul_mod_GPU(inv, div, module);
			int a = mul_mod_GPU(s, matrix[pivot_riga*col+idx], module);
			matrix[idy*col+idx] = sub_mod_GPU(matrix[idy*col+idx], a, module);
		}
	}
}


__global__ void gauss_kernel_celle(int *matrix, int row, int col, int module, int start, int*v, int dim){
		
	int pivot_riga = 0,r = 0,righe_trovate = 0,i,k;
	int s,inv,a;
	int st,flag=0,invarianti=0,flag2=0,tmp;

	if( start == 0 ){
		flag = 1;
	}else{
		st = start;
	}

	for(int pivot_colonna = col-1; pivot_colonna >= 0; pivot_colonna-- ){
		r = righe_trovate;
		while( r < row && matrix[r*col+pivot_colonna] == 0 ){   //m[r][pivot_colonna]
			r++;
			
		}
		// ho trovato la prima riga con elemento non nullo in posizione r e pivot_colonna oppure non esiste nessuna riga con elemento non nullo in posizione pivot_colonna
		
		if( r < row ){ //significa che ho trovato un valore non nullo
			if( r != righe_trovate ){
				swap_rows_GPU(matrix,row,col,righe_trovate,r); //sposto la riga appena trovata nella posizone corretta
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

			inv = invers_GPU(matrix[pivot_riga*col+pivot_colonna],module);		//inverso dell´ elemento in m[r][pivot_colonna]	

			//kernel per riduzione celle
			int block_dim = 16;
			dim3 threads(block_dim, block_dim, 1);
			int block_size = 256;
			int numero_righe = row - righe_trovate;
			int grid_y = numero_righe/block_dim + 1;
			int grid_x = col/block_dim + 1;
			dim3 blocks(grid_x, grid_y,1);

			//printf("Tot celle:%d, Tot thread:%d\n",numero_righe*col,  grid_y*grid_x*block_size);
			
			kernel_riduzione_cella<<<blocks, threads>>>(matrix, row, col, dim, module, inv, pivot_colonna, pivot_riga, righe_trovate);
			cudaDeviceSynchronize();

			//necessario azzerare tutta la colonna (pivot_colonna)
			for(int x = pivot_riga + 1; x < row; x++){
				matrix[x*col+pivot_colonna] = 0;
			}

		}
	}
}


__global__ void gauss_kernel_righe(int *matrix, int row, int col, int module, int start, int*v, int dim){
		
	int pivot_riga = 0,r = 0,righe_trovate = 0,i,k;
	int s,inv,a;
	int st,flag=0,invarianti=0,flag2=0,tmp;


	if( start == 0 ){
		flag = 1;
	}else{
		st = start;
	}

	for(int pivot_colonna = col-1; pivot_colonna >= 0; pivot_colonna-- ){
		r = righe_trovate;
		while( r < row && matrix[r*col+pivot_colonna] == 0 ){   //m[r][pivot_colonna]
			r++;
			
		}
		// ho trovato la prima riga con elemento non nullo in posizione r e pivot_colonna oppure non esiste nessuna riga con elemento non nullo in posizione pivot_colonna
		
		if( r < row ){ //significa che ho trovato un valore non nullo
			if( r != righe_trovate ){
				swap_rows_GPU(matrix,row,col,righe_trovate,r); //sposto la riga appena trovata nella posizone corretta
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

			//printf("Elemento di pivot m[%d][%d]= %d\n", pivot_riga, pivot_colonna, matrix[pivot_riga*col+pivot_colonna]);

			inv = invers_GPU(matrix[pivot_riga*col+pivot_colonna],module);		//inverso dell´ elemento in m[r][pivot_colonna]	

			int block_dim = 1024;
			//kernel per riduzione righe	
			int numero_righe = row - righe_trovate;
			int t = (numero_righe < block_dim ? numero_righe : block_dim);
			//int b = (t < 512 ? 1 : (numero_righe/512) ); 
			int b = 1;			
			if( t == block_dim && numero_righe != block_dim ){
				b = numero_righe / block_dim + 1;
			}

			//printf("Numero di thread lanciati: %d\n", b*t);

			dim3 threads(t);
			dim3 blocks(b);
			kernel_riduzione_riga<<<blocks, threads>>>(matrix, row, col, dim, module, righe_trovate, pivot_colonna, inv, pivot_riga);
			cudaDeviceSynchronize();

		}
	}
}



void gauss_GPU(int *m, int row, int col, int module, int start, int *v){

	int matrix_length = row * col;
	int matrix_length_bytes = matrix_length * sizeof(int);
	
	int *m_d;

	gpuErrchk(cudaMalloc( (void **) &m_d, matrix_length_bytes));
	gpuErrchk(cudaMemcpy(m_d, m, matrix_length_bytes, cudaMemcpyHostToDevice));

	cudaStream_t s;
	cudaStreamCreateWithFlags(&s, cudaStreamNonBlocking);
	cudaDeviceSetLimit(cudaLimitDevRuntimePendingLaunchCount, 30000);
	gauss_kernel_righe<<<1,1,0,s>>>(m_d, row, col, module, start, v, row*col);
	cudaDeviceSynchronize();
	gpuErrchk(cudaMemcpy(m, m_d, matrix_length_bytes, cudaMemcpyDeviceToHost));

	gpuErrchk(cudaFree(m_d));

}

/*  effettua la riduzione di gauss della matrice m in place
    parametro start viene utilizzato per ottimizzare l'agloritmo:
	- 0 effettua la normale riduzione di gauss
	- n (con n > 0) effettua una riduzione ottimizzata saltando le prime n righe nella fase di riduzione delle righe sottostanti.

	Questa ottimizzazione è possibile solo se le prime n righe sono risultato di una riduzione di gauss precedente, quindi risulta inutile ripetere l'operazione di riduzione su queste righe.
	Il salto delle prime n righe avviene solo se nella corrente riduzione non si è effettuato swap di righe. Quando le righe trovate superano start si ritorna alla normale riduzione.
	

*/
void gauss(int *m, int row, int col, int module, int start, int *v){
	
	int pivot_riga = 0,r = 0,righe_trovate = 0,i,k;
	int s,inv,a;
	int st,flag=0,invarianti=0,flag2=0,tmp;

	if( start == 0 ){
		flag = 1;
	}else{
		st = start;
	}

	for(int pivot_colonna = col-1; pivot_colonna >= 0; pivot_colonna-- ){
		r = righe_trovate;
		while( r < row && m[r*col+pivot_colonna] == 0 ){   //m[r][pivot_colonna]
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

			inv = invers(m[pivot_riga*col+pivot_colonna],module);		//inverso dell´ elemento in m[r][pivot_colonna]	

			#pragma omp parallel for private(i,s,k,a)
			for( i = st; i < row; i++ ){
				if( m[i*col+pivot_colonna] != 0 ){		//m[i][pivot_colonna]
					
					s = mul_mod(inv,m[i*col+pivot_colonna],module);						
					for( k = 0; k < pivot_colonna+1; k++ ){
						a = mul_mod(s,m[pivot_riga*col+k],module);
						m[i*col+k] = sub_mod(m[i*col+k],a,module);

					}
				}
			}
		}
	}
}

int null_rows(int *m, int row, int col){
//calcola il numero di righe nulle presenti nella matrice m.

	int i,j,last,null_rows;
	null_rows = 0;
	for(i=0; i<row; i++){
		last = -1;
		for(j=col-1; j>-1; j--){
			if(m[i*col+j] != 0 ){
				last = j;
				break;
			}
		}
		if( last == -1 )
			null_rows++;
	}
	return null_rows;
}

void eliminate_null_rows(int **m, int *row, int col){
//Elimina dalla matrice m le righe nulle.
//N.B. questa procedura elimina le ultime righe nulle della matrice.
//Questa funzione DEVE essere utilizzata dopo la riduzione di Gauss.
//La riduzione di Gauss sposta nelle ultime posizioni tutte le righe nulle.
//Se non si esegue questa funzione dopo Gauss si possono eliminare righe non nulle.	

	int null_row = null_rows(*m,*row,col);
	int new_rows = *row - null_row;
	if(null_row != 0){
		*m = (int *)realloc( *m , (new_rows*col) * sizeof (int));
		*row = new_rows;
	}
}

void print_matrix_degree(int *m_deg, FILE *output_file){
//stampa il vettore dei gradi della matrice.
	int i;
	fprintf(output_file, "Gradi della matrice = {");
	for(i=0; i<max_degree+1; i++)
		if( m_deg[i] != 0 )	fprintf(output_file, " %d ",i);
	fprintf(output_file, "}\n");
}

int target_degree(int *v){
//Controlla se il vettore v rappresenta la condizione di terminazione con gradi completi {1,2,3,...,max_degree}
//Se la condizione è soddisfatta return 0 altrimenti -1

	int i,flag;
	flag = 0;
	for(i=1; i<max_degree+1; i++){
		if( v[i] != 1 ){
			flag = -1;
			break;
		}
	}
	return flag;
}

void execute_standard(int **matrix, int * row, int col, struct map map, int *degree, int **monomi, int numero_variabili, int n_loops, int expansion, FILE *output_file){

	clock_t start, end;
	double elapsed;

	//creo l'array che conterrà i gradi dei vari round
	int **m_deg_array = (int **)malloc(sizeof(int*));
	m_deg_array[0] = (int *)calloc(max_degree+1, sizeof(int));
	int n_round = 0;
	int *m_deg = m_deg_array[0];
	int missing_degree = max_degree;
	fprintf(output_file, "Inizio computazione, metodo standard\n");
	matrix_degree(*matrix, *row, col, m_deg, monomi, numero_variabili);

	int flag, old_v, new_v;
	flag = old_v = new_v = 0;
	old_v = *row;
	//assumo espansione normale
	int expansion_degree = max_degree;
	int st = 0;
	
	while( flag != 1 ){
		n_round++;

		fprintf(output_file, "\n -Eseguo moltiplicazione, ");
		fflush(stdout);

		start = clock();	

		//find missing degree to multiply matrix
			for(int i=max_degree; i>0; i--){
				if( m_deg[i] == 0 ){
					missing_degree = i;
					break;
				}
			}

		moltiplica_matrice(matrix, row, col, map, degree, monomi, numero_variabili, missing_degree);

		end = clock();
		elapsed =  ((double)(end - start)) / CLOCKS_PER_SEC;
		fprintf(output_file, "numero righe: %d     (%f sec)", *row, elapsed);

		fprintf(output_file, "\n -Eseguo Gauss, ");
		fflush(stdout);
		start = clock();

		//applico la riduzione di Gauss
		//gauss(*matrix, *row, col, module, st, NULL);
		gauss_GPU(*matrix, *row, col, module, st, NULL);
		//elimino le righe nulle della matrice
		eliminate_null_rows(matrix, row, col);
		
		//aggiungo all'array i gradi dell'attuale round
		//n_round+1 perchè salvo anche i gradi prima di inziare i round
		m_deg_array = (int **)realloc(m_deg_array, sizeof(int*)*(n_round+1));
		m_deg_array[n_round] = (int *)calloc(max_degree+1, sizeof(int));
		m_deg = m_deg_array[n_round];
		
		end = clock();
		elapsed =  ((double)(end - start)) / CLOCKS_PER_SEC;
		fprintf(output_file, "numero righe: %d               (%f sec)\n", *row, elapsed);

  		matrix_degree(*matrix,*row, col, m_deg, monomi, numero_variabili);
		print_matrix_degree(m_deg, output_file);

		new_v = *row;
		st = new_v;


		if( target_degree(m_deg) == 0 )
			flag = 1;
		else{
			old_v = new_v;
			}

		//flag = 1;
	/*	
		if(n_round == 2){
			flag=1;
			//print_matrix(*matrix, *row, col, output_file);
		}
	*/	
		for(int i=max_degree; i>0; i--){
			if( m_deg[i] == 0 ){
				expansion_degree = i;
				break;
			}
		}


	}
	for (int i = 0; i < n_round+1; i++)
		free(m_deg_array[i]);	
	free(m_deg_array);	
}


void print_incognite(int *m, int row, int col, int num_var, int **vet, FILE *output_file){

	int grado,last;

	for(int r = row - (num_var+1); r<row; r++){

		//cerca la posizione dell'ulitmo elemento non nullo della riga r
		for( int i=col-1; i>=0; i-- ){
			if( m[r*col+i] != 0 ){ //m[r][i] != 0
				last = i;
				break;
			}
		}
		//calcola il grado della riga r
		grado = grado_monomio(last,vet,num_var);
		//se il grado della riga r è 1 allora stampa tutta la riga della matrice
		if( grado == 1 ){
			for( int j=0; j<last+1; j++ ){
				fprintf(output_file, "%d ", m[r*col+j]); //m[r][j]
			}
			fprintf(output_file, "\n\n");
		}
	}
	fprintf(output_file, "\n");	
}

int main (int argc, char *argv[]){

	FILE *input_file = NULL, *output_file = NULL;
	//inizializzo flag a false
	
	for (int parsed = 1; parsed < argc; parsed++) {
		if (parsed < argc && !strcmp(argv[parsed], "--input")) {
			parsed++;
			input_file = fopen(argv[parsed], "r");
			if (!input_file) {
				perror("Errore nell'apertura del file di input");
				return (-1);
			}
		}
		else if (parsed < argc && !strcmp(argv[parsed], "--output")) {
			parsed++;
			output_file = fopen(argv[parsed], "w");
			if (!output_file) {
				perror("Errore nell'apertura del file di output");
				return (-1);
			}
		}
	}


	if (!input_file)
		input_file = stdin;
	if (!output_file)
		output_file = stdout;


	int row, col, numero_variabili, tipo_ordinamento;
	int *matrix;

	char *variabili;
	row = col = numero_variabili = 0;
	int (*ord) (void*, const void *, const void *);
	int *d_row, **map;
	struct map smap;

	clock_t start, end;
	double elapsed = 0.0;
	start = clock();

	//alloca la matrice principale, legge da input: il modulo,massimo grado e numero variabili
	allocation(&matrix, &row, &col, &numero_variabili, &variabili, &tipo_ordinamento, &module, &max_degree, input_file);
	
	int matrix_lentght = row * col; //numero di elementi della matrice

	if( order(&ord, tipo_ordinamento) != 0 ){
		fprintf(stderr, "Ordinamento insesistente!!!\n\nTERMINAZIONE PROGRAMMA");
		return 0;
	}

	int * degree = (int *)calloc(max_degree+1, sizeof(int));
	int numero_monomi = col;
	int **monomi;

	//crea il vettore con tutti i possibili monomi avendo num_var varaibili e max_degree come massimo grado
	monomi = monomial_computation(numero_variabili, max_degree, numero_monomi);

	//ordina il vettore dei monomi secondo un determinato ordinamento, ordinamento intercambiabile
	qsort_s(monomi, numero_monomi, sizeof(int*), ord, &numero_variabili);
		
	//inizializzazione matrice (lettura dati input)
	if (init_matrix(matrix, row, col, monomi, variabili, numero_variabili, ord, input_file) == -1) {
		fprintf(stderr, "Errore di input !!!\n\nTERMINAZIONE PROGRAMMA"); //se l'input è in formato scorrettro abort del programma
		return 0;
	}
	
	end = clock();
	elapsed =  ((double)(end - start)) / CLOCKS_PER_SEC;
	fprintf(output_file, "\nInizializzazione in %f sec\n", elapsed);	

	start = clock();
	
	setup_struct_map(&smap, monomi, numero_monomi, numero_variabili, max_degree, ord);
	
	end = clock();
	elapsed =  ((double)(end - start)) / CLOCKS_PER_SEC;
	fprintf(output_file, "\nMappa creata in %f sec,   %d x %d \n\n", elapsed, col, col);

	//RISOLUZIONE PROBLEMA
	start = clock();	

	//inizializzazione vettore dei gradi dei polinomi
	init_degree_vector(degree, numero_variabili);

	int n_loops = 30, expansion = 1;
	//eseguo moltiplicazione e riduzione di Gauss finche non trovo soluzione
	execute_standard(&matrix, &row, col, smap, degree, monomi, numero_variabili, n_loops, expansion, output_file);

	end = clock();
	elapsed =  ((double)(end - start)) / CLOCKS_PER_SEC;
	fprintf(output_file, "\nTarget raggiunto, soluzione trovata in %f sec\n\n", elapsed);

	//print_matrix(matrix, row, col, output_file);
	print_incognite(matrix, row, col, numero_variabili, monomi, output_file);
	for(int i=0; i<row*col; i++){
		if(matrix[i] > module){
			printf("OVERFLOW\n");
		}
	}

	free(matrix);
	free(degree);
	cudaDeviceReset();

	return 0;
}

