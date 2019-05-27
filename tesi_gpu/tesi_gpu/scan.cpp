#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include "linalg.h"
#include "matrix.h"


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
	*variabili = (char *) malloc(sizeof(char));
	c = fgetc(input_file);
	while (c != '\n') {
		(*variabili)[i] = c;
		i++;
		(*numero_variabili)++;
		*variabili = (char *) realloc(*variabili, (i + 1) * sizeof(char));
		c = fgetc(input_file);
	}

	*col = monomial_combinations(*numero_variabili, *max_degree);
	/*
	*matrix = malloc((*row) * sizeof (int *) );            // allocazione della matrice dei coefficienti
	if( *matrix != NULL )
	for (int i=0; i<(*row); i++)
	(*matrix)[i] = calloc((*col) , sizeof (int) );
	*/
	*matrix = (int *)calloc((*row) * (*col), sizeof(int));

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

		char *coefficient = (char *) malloc(sizeof(char));
		while (index < len && isdigit(mon[index])) {
			coefficient = (char *) realloc(coefficient, (index + 1) * sizeof(char));
			coefficient[index] = mon[index];
			index++;
		}
		//aggiungo il carattere di temrinazione
		coefficient = (char *) realloc(coefficient, (index + 1) * sizeof(char));
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
			exponent = (char *) malloc(sizeof(char));
			while (index < len && isdigit(mon[index])) {
				exponent = (char *) realloc(exponent, (exponent_index + 1) * sizeof(char));
				exponent[exponent_index] = mon[index];
				exponent_index++;
				index++;
			}
			//metto il carattere di terminazoine stringa
			exponent = (char *) realloc(exponent, (exponent_index + 1) * sizeof(char));
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







int order(int(**ord) (void *, const void *, const void*), int n) {
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


