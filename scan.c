#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include "scan.h"
#include "linalg.h"

void allocation(long long ***m, int *row, int *col, int *num_var, char **v, int *n, long long *module, int *max_degree, FILE *input_file){
/*
Legge da input le seguenti informazioni:
	- modulo dei coefficienti
	- grado massimo
	- numero dei polinomi di partenza
	- tipo di ordinamento
	- variabili utilizzate nei polinomi


con queste informazioni alloca la matrice principale (matrice che conterrà i polinomi) e stabilisce il numero di variabili utilizzate.
*/
	fscanf(input_file, "%lli",module); //leggo il modulo
	fgetc(input_file);
	fscanf(input_file, "%d",max_degree); //leggo il grado massimo
	fgetc(input_file);
	fscanf(input_file, "%d",row);  //leggo numero dei polinomi di partenza
	fgetc(input_file);
	fscanf(input_file, "%d",n);  //leggo tipo di ordinamento
	fgetc(input_file);

	int i,j,k,pos_pol,num_pol;
	char c;

	i=0;
	pos_pol = 0;
	*v = malloc(sizeof(char));
	c = fgetc(input_file);
	while( c != '\n' ){
		(*v)[i] = c;
		i++;
		(*num_var)++;
		*v = realloc(*v, (i+1)*sizeof(char) );
		c = fgetc(input_file);
	}

	*col = monomial_combinations(*num_var, *max_degree);

	*m = malloc((*row) * sizeof (long long *) );            // allocazione della matrice dei coefficienti
	if( *m != NULL )
		for (int i=0; i<(*row); i++)
			(*m)[i] = calloc((*col) , sizeof (long long) );	


}


int parse(int num_var, char *vet, long long **m, int row, int **vet_grd, int len, long long module, int (*ord) (const void *, const void *, void*), FILE *input_file){
/*
Esegue la lettura (parse) dei polinomi di partenza nel seguente modo.
Si legge un monomio alla volta. 
Il monomio viene scomposta dalla funzione parse_mon.
Si inserisce il coefficiente del monomio nella matrice principale (matrice dei coefficienti) nella posizione corretta.
La posizione corretta è indicata da vet_grd.
Si leggono tutti i monomi di tutti i polinomi di partenza.
In caso di errore di formato nell'input la funzione si interrompe restituendo segnale di errore -1.
*/
	int pos_pol = 0,i,col;
	char c,* mon;
	long long cof = 0;
	c = fgetc(input_file);

	int *grade;

	//finchè non termino il file o non ho terminato il numero di polinomi dichiarati
	while( c != EOF && pos_pol < row ){
		mon = malloc( sizeof(char) );
		grade = calloc(num_var,sizeof(int));
		i = 0;	
		while( c != '+' && c != EOF && c != '\n'){
			mon = realloc(mon, (i+1)* sizeof(char));
			mon[i] = c;
			i++;
			c = fgetc(input_file);
		}
		//se non ho salvato niente in mon (i = 0) non faccio il parsing
		if(i != 0 && parse_mon(mon,i,&cof,num_var,vet,grade,module) == -1 ){
			return -1;
		}
		//inserire monomio in posizione corretta
		col = (int **)(bsearch_r((void *) &grade, (void *) vet_grd, len, (sizeof(int*)), ord, &num_var)) - vet_grd;
		m[pos_pol][col] = cof;
		if(c=='\n'){
			pos_pol++;
		}
		free(mon);
		free(grade);
		c = fgetc(input_file);
	}
	return 0;
}


/* mon: stringa che rappresenta un monomio (non c'è carattere terminazione stringa)
 * len: numero di caratteri in mon
 * val: variabile in cui salvare il coefficiente del monomio
 * num_var: numero di variabili nel sistema
 * vet: vettore di caratteri in cui ogni carattere è una variabile (letto precedentemente da input)
 * grade: vettore in cui salvare i gradi delle variabili secondo l'ordine di vet
 * module: campo su cui è rappresentato il sistema 
 */
int parse_mon(char * mon, int len,long long * val, int num_var, char *vet, int *grade, long long module){

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

		char *coefficient = malloc(sizeof(char));
		while (index < len && isdigit(mon[index])) {
			coefficient = realloc(coefficient, (index+1)*sizeof(char));
			coefficient[index] = mon[index];
			index++;
		}
		//aggiungo il carattere di temrinazione
		coefficient = realloc(coefficient, (index+1)*sizeof(char));
		coefficient[index] = '\0';
		//traduco il coefficiente in valore numerico e calcolo il modulo
		*val = mod(atoll(coefficient),module);  
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
			exponent = malloc(sizeof(char));
			while (index < len && isdigit(mon[index])) {
				exponent = realloc(exponent, (exponent_index+1)*sizeof(char));
				exponent[exponent_index] = mon[index];
				exponent_index++;
				index++;
			}
			//metto il carattere di terminazoine stringa
			exponent = realloc(exponent, (exponent_index+1)*sizeof(char));
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




int order(int (**ord) (const void *, const void *, void*), int n){
//inizializza il puntatore ord alla funzione di ordinamento adeguata. Il numero n indica quale funzione scegliere.

	switch(n){

		case 0:
			*ord = grevlex_comparison;
			return 0;
			break;

		default:
			return -1;
			break;	

	}
}


