/*
* 
* Software necessario per l'esecuzione:
* 
*    NASM (www.nasm.us)
*    GCC (gcc.gnu.org)
* 
* entrambi sono disponibili come pacchetti software 
* installabili mediante il packaging tool del sistema 
* operativo; per esempio, su Ubuntu, mediante i comandi:
* 
*    sudo apt-get install nasm
*    sudo apt-get install gcc
* 
* potrebbe essere necessario installare le seguenti librerie:
* 
*    sudo apt-get install lib64gcc-4.8-dev (o altra versione)
*    sudo apt-get install libc6-dev-i386
* 
* Per generare il file eseguibile:
* 
* nasm -f elf64 fss64.nasm && gcc -m64 -msse -O0 -no-pie sseutils64.o fss64.o fss64c.c -o fss64c -lm && ./fss64c $pars
* 
* oppure
* 
* ./runfss64
* 
*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <libgen.h>
#include <xmmintrin.h>
#include <stdbool.h>

#include <omp.h>

#define	type		double
#define	MATRIX		type*
#define	VECTOR		type*

#define UNROLL_FACTOR 4

typedef struct {
	MATRIX ds; 		// dataset
	VECTOR labels; 	// etichette
	int* out;		// vettore contenente risultato dim=k
	type sc;		// score dell'insieme di features risultato
	int k;			// numero di features da estrarre
	int N;			// numero di righe del dataset
	int d;			// numero di colonne/feature del dataset
	int display;
	int silent;
} params;

/*
* 
*	Le funzioni sono state scritte assumento che le matrici siano memorizzate 
* 	mediante un array (float*), in modo da occupare un unico blocco
* 	di memoria, ma a scelta del candidato possono essere 
* 	memorizzate mediante array di array (float**).
* 
* 	L'assunzione corrente è che le matrici siano lette dal file in row-major order,
*   ma poi trasposte per averle in column-major order.
*
* 
*/

void* get_block(int size, int elements) {
    return _mm_malloc(elements*size,32);
}

void free_block(void* p) { 
	_mm_free(p);
}

MATRIX alloc_matrix(int rows, int cols) {
	return (MATRIX) get_block(sizeof(type),rows*cols);
}

int* alloc_int_matrix(int rows, int cols) {
	return (int*) get_block(sizeof(int),rows*cols);
}

void dealloc_matrix(void* mat) {
	free_block(mat);
}

/*
* 
* 	load_data
* 	=========
* 
*	Legge da file una matrice di N righe
* 	e M colonne e la memorizza in un array lineare in row-major order
* 
* 	Codifica del file:
* 	primi 4 byte: numero di colonne (N) --> numero intero
* 	successivi 4 byte: numero di righe (M) --> numero intero
* 	successivi N*M*8 byte: matrix data in row-major order --> numeri floating-point a precisione doppia
* 
* 
*/
MATRIX load_data(char* filename, int *n, int *k) {
	FILE* fp;
	int rows, cols, status, i;
	
	fp = fopen(filename, "rb");
	
	if (fp == NULL){
		printf("'%s': bad data file name!\n", filename);
		exit(0);
	}
	
	status = fread(&cols, sizeof(int), 1, fp);
	status = fread(&rows, sizeof(int), 1, fp);
	
	MATRIX data = alloc_matrix(rows,cols);
	status = fread(data, sizeof(type), rows*cols, fp);
	fclose(fp);
	
	*n = rows;
	*k = cols;
	
	return data;
}

/*
* 	save_data
* 	=========
* 
*	Salva su file un array lineare in row-major order
*	come matrice di N righe e M colonne
* 
* 	Codifica del file:
* 	primi 4 byte: numero di colonne (N) --> numero intero a 32 bit
* 	successivi 4 byte: numero di righe (M) --> numero intero a 32 bit
* 	successivi N*M*8 byte: matrix data in row-major order --> numeri interi o floating-point a precisione doppia
*/
void save_data(char* filename, void* X, int n, int k) {
	FILE* fp;
	int i;
	fp = fopen(filename, "wb");
	if(X != NULL){
		fwrite(&k, 4, 1, fp);
		fwrite(&n, 4, 1, fp);
		for (i = 0; i < n; i++) {
			fwrite(X, sizeof(type), k, fp);
			//printf("%i %i\n", ((int*)X)[0], ((int*)X)[1]);
			X += sizeof(type)*k;
		}
	}
	else{
		int x = 0;
		fwrite(&x, 4, 1, fp);
		fwrite(&x, 4, 1, fp);
	}
	fclose(fp);
}

/*
* 	save_out
* 	=========
* 
*	Salva su file un array lineare composto da k+1 elementi.
* 
* 	Codifica del file:
* 	primi 4 byte: contenenti il numero di elementi (k+1)		--> numero intero a 32 bit
* 	successivi 4 byte: numero di righe (1) 						--> numero intero a 32 bit
* 	successivi byte: elementi del vettore 		--> 1 numero floating-point a precisione doppia e k interi
*/
void save_out(char* filename, type sc, int* X, int k) {
	FILE* fp;
	int i;
	int n = 1;
	k++;
	fp = fopen(filename, "wb");
	if(X != NULL){
		fwrite(&n, 4, 1, fp);
		fwrite(&k, 4, 1, fp);
		fwrite(&sc, sizeof(type), 1, fp);
		fwrite(X, sizeof(int), k, fp);
		//printf("%i %i\n", ((int*)X)[0], ((int*)X)[1]);
	}
	fclose(fp);
}

// PROCEDURE ASSEMBLY
extern void correlation_cf(type* feature, type* labels, int n, type* r_cf);
extern void correlation_ff(type* x, type* y, int n, type* r_ff);
extern void transpose(MATRIX A, MATRIX B, int ncols, int nrows);

//Funzione principale di Correlation Features Selection (CFS) 
void cfs(params* input) {
    
    //Dichiarazione variabili                                                                           
    int selectedFeatures[input->k]; //Array per tracciare le caratteristiche selezionate 
    bool featureSelected[input->d]; //Array booleano per marcare se una caratteristica è stata selezionata
    memset(featureSelected, 0, sizeof(featureSelected)); //Inizializzazione di featureSelected a false
    
    type cf_vector[input->d]; //Array per la correlazione class-feature
    type **ff_matrix = (type **)malloc(input->d * sizeof(type *)); //Matrice per la correlazione feature-feature
    #pragma omp parallel for
    for(int i = 0; i < input->d; i++) {
        ff_matrix[i] = (type *)malloc(input->d * sizeof(type));
    }
    
    type r_cf = 0.0; //Correlazione class-feature
    type r_cf_sum = 0.0; //Somma delle correlazioni class-feature
    type r_cf_den = 0.0; //Denominatore per il calcolo del merito
    type r_ff = 0.0; //Correlazione feature-feature
    type r_ff_sum = 0.0; //Somma delle correlazioni feature-feature
    type r_ff_den = 0.0; //Denominatore per il calcolo del merito
    
    #pragma omp parallel for private(r_cf)
    for (int i = 0; i < input->d; i++) {
        correlation_cf(input->ds + i * input->N, input->labels, input->N, &r_cf);
        cf_vector[i] = r_cf;
 
        // Unrolling del loop interno
        int j;
        for (j = 0; j <= input->d - UNROLL_FACTOR; j += UNROLL_FACTOR) {
            ff_matrix[i][j] = -2.0;
            ff_matrix[i][j + 1] = -2.0;
            ff_matrix[i][j + 2] = -2.0;
            ff_matrix[i][j + 3] = -2.0;
        }
 
        // Gestione dei casi in cui input->d non è un multiplo di UNROLL_FACTOR
        for (; j < input->d; j++) {
            ff_matrix[i][j] = -2.0;
        }
    }
    
    //Ciclo di selezione delle caratteristiche
    for(int s = 0; s < input->k; s++) { 
        type maxMerit = -1.0; //Variabile per tenere traccia del merito massimo
        int maxFeature = -1; //Variabile per tenere traccia della miglior caratteristica
        type r_cf_sum_tmp = 0.0; 
        type r_ff_sum_tmp = 0.0; 
        r_cf_den += 1.0; //Aggiornamento del denominatore class-feature
        r_ff_den += (type) s; //Aggiornamento del denominatore feature-feature
        type r_cf_sum_max = r_cf_sum; 
        type r_ff_sum_max = r_ff_sum;
        
        //Ciclo interno per valutare ogni caratteristica non selezionata  
        for(int i = 0; i < input->d; i++) { 
            if(!featureSelected[i]) {
                r_cf_sum_tmp = r_cf_sum; 
                r_ff_sum_tmp = r_ff_sum;  
                r_cf = cf_vector[i]; 
                r_cf_sum_tmp += fabs(r_cf);
                
                //Calcolo delle correlazioni feature-feature per la caratteristica corrente               
                for(int j = 0; j < s; j++) { 
                    if(ff_matrix[i][j] != -2.0) {
                        r_ff_sum_tmp += ff_matrix[i][j];
                    }
                    else {
                        correlation_ff(input->ds + i * input->N, input->ds + selectedFeatures[j] * input->N, input->N, &r_ff);
                        type t = fabs(r_ff);
                        r_ff_sum_tmp += t;
                        ff_matrix[i][j] = t;
                    }
                }
                
                //Calcolo del merito della caratteristica corrente 
                type K = (type) s + 1.0; 
                type merit = (s == 0) ? (K  * fabs(r_cf))/sqrt(K) : 
                (K * (r_cf_sum_tmp / r_cf_den)) / sqrt(K + (K * (K - 1.0) * (r_ff_sum_tmp / r_ff_den)));
                
                //Aggiornamento del massimo merito e della caratteristica corrispondente   
                if (merit > maxMerit) { 
                    maxMerit = merit; 
                    maxFeature = i;
                    r_cf_sum_max = r_cf_sum_tmp;
                    r_ff_sum_max = r_ff_sum_tmp;
                } 
            } 
        }
        
        //Aggiornamento delle variabili dopo la selezione della caratteristica 
        r_cf_sum = r_cf_sum_max;
        r_ff_sum = r_ff_sum_max;
        featureSelected[maxFeature] = true; 
        selectedFeatures[s] = maxFeature; 
        input->sc = maxMerit;            
        input->out[s] = maxFeature; 
    }
    
    #pragma omp parallel for
    for (int i = 0; i < input->d; i++) { //Rilascio memoria
        free(ff_matrix[i]);
    }
    free(ff_matrix);

}

int main(int argc, char** argv) {

	char fname[256];
	char* dsfilename = NULL;
	char* labelsfilename = NULL;
	clock_t t;
	float time;
	
	//
	// Imposta i valori di default dei parametri
	//

	params* input = malloc(sizeof(params));

	input->ds = NULL;
	input->labels = NULL;
	input->k = -1;
	input->sc = -1;

	input->silent = 0;
	input->display = 0;

	//printf("%i\n", sizeof(int));

	//
	// Visualizza la sintassi del passaggio dei parametri da riga comandi
	//

	if(argc <= 1){
		printf("%s -ds <DS> -labels <LABELS> -k <K> [-s] [-d]\n", argv[0]);
		printf("\nParameters:\n");
		printf("\tDS: il nome del file ds2 contenente il dataset\n");
		printf("\tLABELS: il nome del file ds2 contenente le etichette\n");
		printf("\tk: numero di features da estrarre\n");
		printf("\nOptions:\n");
		printf("\t-s: modo silenzioso, nessuna stampa, default 0 - false\n");
		printf("\t-d: stampa a video i risultati, default 0 - false\n");
		exit(0);
	}

	//
	// Legge i valori dei parametri da riga comandi
	//

	int par = 1;
	while (par < argc) {
		if (strcmp(argv[par],"-s") == 0) {
			input->silent = 1;
			par++;
		} else if (strcmp(argv[par],"-d") == 0) {
			input->display = 1;
			par++;
		} else if (strcmp(argv[par],"-ds") == 0) {
			par++;
			if (par >= argc) {
				printf("Missing dataset file name!\n");
				exit(1);
			}
			dsfilename = argv[par];
			par++;
		} else if (strcmp(argv[par],"-labels") == 0) {
			par++;
			if (par >= argc) {
				printf("Missing labels file name!\n");
				exit(1);
			}
			labelsfilename = argv[par];
			par++;
		} else if (strcmp(argv[par],"-k") == 0) {
			par++;
			if (par >= argc) {
				printf("Missing k value!\n");
				exit(1);
			}
			input->k = atoi(argv[par]);
			par++;
		} else{
			printf("WARNING: unrecognized parameter '%s'!\n",argv[par]);
			par++;
		}
	}

	//
	// Legge i dati e verifica la correttezza dei parametri
	//

	if(dsfilename == NULL || strlen(dsfilename) == 0){
		printf("Missing ds file name!\n");
		exit(1);
	}

	if(labelsfilename == NULL || strlen(labelsfilename) == 0){
		printf("Missing labels file name!\n");
		exit(1);
	}


	input->ds = load_data(dsfilename, &input->N, &input->d);
	MATRIX dataset_cmo = alloc_matrix(input->d, input->N);
	transpose(input->ds, dataset_cmo, input->d, input->N);
	dealloc_matrix(input->ds);
	input->ds = dataset_cmo;

	int nl, dl;
	input->labels = load_data(labelsfilename, &nl, &dl);
	
	if(nl != input->N || dl != 1){
		printf("Invalid size of labels file, should be %ix1!\n", input->N);
		exit(1);
	} 

	if(input->k <= 0){
		printf("Invalid value of k parameter!\n");
		exit(1);
	}

	input->out = alloc_int_matrix(input->k, 1);

	//
	// Visualizza il valore dei parametri
	//

	if(!input->silent){
		printf("Dataset file name: '%s'\n", dsfilename);
		printf("Labels file name: '%s'\n", labelsfilename);
		printf("Dataset row number: %d\n", input->N);
		printf("Dataset column number: %d\n", input->d);
		printf("Number of features to extract: %d\n", input->k);
	}


	//
	// Correlation Features Selection
	//
	t = clock();
	cfs(input);
	t = clock() - t;
	time = ((float)t)/CLOCKS_PER_SEC;

	if(!input->silent)
		printf("CFS time = %.3f secs\n", time);
	else
		printf("%.3f\n", time);

	//
	// Salva il risultato
	//
	sprintf(fname, "out64_%d_%d_%d.ds2", input->N, input->d, input->k);
	save_out(fname, input->sc, input->out, input->k);
	if(input->display){
		if(input->out == NULL)
			printf("out: NULL\n");
		else{
			int i,j;
			printf("sc: %lf, out: [", input->sc);
			for(i=0; i<input->k; i++){
				printf("%i,", input->out[i]);
			}
			printf("]\n");
		}
	}

	if(!input->silent)
		printf("\nDone.\n");

	dealloc_matrix(input->ds);
	dealloc_matrix(input->labels);
	dealloc_matrix(input->out);
	free(input);
	
	return 0;
}
