/*
 * corpus.h
 *
 *	Note: word is an instance of term
 *
 */

#ifndef CORPUS_H_
#define CORPUS_H_

#include <gsl/gsl_permutation.h>

typedef struct para {
	int total;
	int* words;
	int* words_local;
} para;

typedef struct doc {
	int total;
	int nparas;
	para* paras;
} doc;

typedef struct Corpus {
	int I;
	int totalw;
	int* total;
	int* ndocs;
	int* max_length; // the maximum length of paragraphs
	int* W_i;
	int** l2g_vocabulary; // local to global vocabulary index
	doc** docs;
} Corpus;

typedef struct id_2_word {
	int id;
	char word_str[100];
} id_2_word;

typedef struct vocabulary {
	int size;
	id_2_word* word_map;
} vocabulary;

/*
 * functions
 */
void read_data(char* dir, char** file, int I, Corpus* c_tr, Corpus* c_te, double percent);
void write_corpus(Corpus*, char**);

void save_local_vocabulary(Corpus*, vocabulary*, char*);

void split_corpus(Corpus* s, Corpus* c1, Corpus* c2);

vocabulary* read_vocabulary(char*);
void write_vocabulary(vocabulary*, char*);

void free_corpus(Corpus*);
void free_vocabulary(vocabulary*);

#endif /* CORPUS_H_ */
