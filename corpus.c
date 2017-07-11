/*
 * corpus.c
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <sys/stat.h>

#include "corpus.h"

int v_size;

void read_data(char* dir, char** file, int I, Corpus* c_tr, Corpus* c_te, double percent)
{
	FILE* fileptr;
	Corpus* c;
	char fileName[BUFSIZ];
	int i, j, l, ndocs, nparas, nwords;
	int word;
	int wtotal = 0;
	int ptotal = 0;
	int max_length = 0;
	int tmp;
	c = (Corpus*) malloc(sizeof(Corpus));
	c->I = I;
	c->max_length = (int*)malloc(sizeof(int) * I);
	c->ndocs = (int*)malloc(sizeof(int) * I);
	c->total = (int*)malloc(sizeof(int) * I);
	c->W_i = (int*)malloc(sizeof(int) * I);
	c->docs = (doc**)malloc(sizeof(doc*) * I);
	c->l2g_vocabulary = (int**)malloc(sizeof(int*) * I);

	c->totalw = 0;

	for(i = 0; i < I; ++i){
		printf("reading data from %s/%s\n", dir, file[i]);
		sprintf(fileName, "%s/%s", dir, file[i]);
		fileptr = fopen(fileName, "r");
		if(!fileptr){
			printf("Cannot open file %s/%s\n", dir, file[i]);
			exit(0);
		}
		ndocs = 0;
		ptotal = 0;
		wtotal = 0;
		max_length = 0;
		c->docs[i] = (doc*) malloc(sizeof(doc) * 1);
		while (1){
			tmp = fscanf(fileptr, "doc:%d\n", &nparas);
			assert(tmp != 0);
			if(tmp == EOF)
				break;
			c->docs[i] = (doc*) realloc(c->docs[i], sizeof(doc) * (ndocs + 1));
			doc* d = &(c->docs[i][ndocs]);
			d->total = 0;
			d->nparas = 1;
			d->paras = (para*) malloc(sizeof(para) * d->nparas);

			para* p = &(d->paras[0]);
			p->words = (int*)malloc(sizeof(int) * 1);
			p->words_local = (int*)malloc(sizeof(int) * 1);
			p->total = 0;
			for (j = 0; j < nparas; j++){
				tmp = fscanf(fileptr, "para:%d\n", &nwords);
				assert(tmp != 0);
				p->total += nwords;

				if(nwords > max_length){
					max_length = nwords;
				}

				p->words = (int*)realloc(p->words, sizeof(int) * p->total);
				p->words_local = (int*)realloc(p->words_local, sizeof(int) * p->total);
				for (l = 0; l < nwords; l++){
					tmp = fscanf(fileptr, "word:%d\n", &word);
					assert(tmp != 0);
					if(word >= v_size){
						printf("warning: word index exceeds vocabulary size, force it to the max index allowed!!!\n");
						word = v_size - 1;
					}
					p->words[l + p->total - nwords] = word;
#ifdef FL
					p->words_local[l + p->total - nwords] = word;
#endif
				}
				d->total += nwords;
				ptotal += 1;
			}
			wtotal += d->total;
			ndocs++;
		}
		fclose(fileptr);
		c->ndocs[i] = ndocs;
		c->total[i] = wtotal;
		c->max_length[i] = max_length;
		c->totalw += wtotal;

		int id, nW = 0;
#ifdef FL
		nW = v_size;
#else
		int tag, it;
		int *wordInd = (int*)malloc(sizeof(int) * wtotal);
		for(id = 0; id < ndocs; ++id){
			for(it = 0; it < c->docs[i][id].paras[0].total; ++it){
				int wordd = c->docs[i][id].paras[0].words[it];
				tag = 0;
				for(j = 0; j < nW; ++j){
					if(wordd == wordInd[j]){
						tag = 1;
						break;
					}
				}
				if(tag == 0){
					wordInd[nW] = wordd;
					++nW;
				}
				c->docs[i][id].paras[0].words_local[it] = j;
			}
		}
#endif
		c->W_i[i] = nW;
		c->l2g_vocabulary[i] = (int*)malloc(sizeof(int) * nW);
		for(id = 0; id < nW; ++id){
#ifdef FL
			c->l2g_vocabulary[i][id] = id;
#else
			c->l2g_vocabulary[i][id] = wordInd[id];
#endif
		}
#ifndef FL
		free(wordInd);
#endif
		printf("\n>>>>>>statistic for loaded corpus <<<<<<\n");
		printf("number of documents      : %d\n", ndocs);
		printf("paragraphs total         : %d\n", ptotal);
		printf("words total              : %d\n", wtotal);
		printf("max length               : %d\n", max_length);
		printf("local vocabulary size    : %d\n", nW);
		printf(">>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<\n");
	}

	int num_tr, num_te, i1;
	c_tr->I = c->I;
	c_tr->totalw = 0;
	c_tr->total = (int*)malloc(sizeof(int) * c->I);
	c_tr->ndocs = (int*)malloc(sizeof(int) * c->I);
	c_tr->max_length = (int*)malloc(sizeof(int) * c->I);
	c_tr->W_i = (int*)malloc(sizeof(int) * c->I);
	c_tr->l2g_vocabulary = (int**)malloc(sizeof(int*) * c->I);
	c_tr->docs = (doc**)malloc(sizeof(doc*) * c->I);
	if(percent != 1){
		c_te->I = c->I;
		c_te->totalw = 0;
		c_te->total = (int*)malloc(sizeof(int) * c->I);
		c_te->ndocs = (int*)malloc(sizeof(int) * c->I);
		c_te->max_length = (int*)malloc(sizeof(int) * c->I);
		c_te->W_i = (int*)malloc(sizeof(int) * c->I);
		c_te->l2g_vocabulary = (int**)malloc(sizeof(int*) * c->I);
		c_te->docs = (doc**)malloc(sizeof(doc*) * c->I);
	}

	for(i = 0; i < c->I; i++){
		assert(c->ndocs[i] > 1);
		num_tr = (int)(c->ndocs[i] * percent);
		if(percent != 1 && num_tr == c->ndocs[i]){
			num_tr--;
		}
		num_te = c->ndocs[i] - num_tr;
		c_tr->docs[i] = (doc*)malloc(sizeof(doc)*num_tr);
		c_tr->ndocs[i] = num_tr;
		c_tr->total[i] = 0;
		c_tr->max_length[i] = 0;
		c_tr->W_i[i] = c->W_i[i];
		c_tr->l2g_vocabulary[i] = (int*)malloc(sizeof(int) * c->W_i[i]);
		for(i1 = 0; i1 < c->W_i[i]; i1++){
			c_tr->l2g_vocabulary[i][i1] = c->l2g_vocabulary[i][i1];
		}
		if(percent != 1){
			c_te->docs[i] = (doc*)malloc(sizeof(doc)*num_te);
			c_te->ndocs[i] = num_te;
			c_te->total[i] = 0;
			c_te->max_length[i] = 0;
			c_te->W_i[i] = c->W_i[i];
			c_te->l2g_vocabulary[i] = (int*)malloc(sizeof(int) * c->W_i[i]);
			for(i1 = 0; i1 < c->W_i[i]; i1++){
				c_te->l2g_vocabulary[i][i1] = c->l2g_vocabulary[i][i1];
			}
		}

		for(i1 = 0; i1 < c_tr->ndocs[i]; i1++){
			c_tr->docs[i][i1].paras = (para*)malloc(sizeof(para)*c->docs[i][i1].nparas);
			c_tr->docs[i][i1].nparas = c->docs[i][i1].nparas;
			c_tr->docs[i][i1].total = c->docs[i][i1].total; //i1 different in c_te
			c_tr->totalw += c_tr->docs[i][i1].total;

			for(j = 0; j < c_tr->docs[i][i1].nparas; j++){
				c_tr->docs[i][i1].paras[j].words = (int*)malloc(sizeof(int) * c->docs[i][i1].paras[j].total); //i1 different in c_te
				c_tr->docs[i][i1].paras[j].words_local = (int*)malloc(sizeof(int) * c->docs[i][i1].paras[j].total); //i1 different in c_te
				c_tr->docs[i][i1].paras[j].total = c->docs[i][i1].paras[j].total;
				c_tr->total[i] += c->docs[i][i1].paras[j].total;
				for(l = 0; l < c->docs[i][i1].paras[j].total; l++){ //i1 different in c_te
					c_tr->docs[i][i1].paras[j].words[l] = c->docs[i][i1].paras[j].words[l];
					c_tr->docs[i][i1].paras[j].words_local[l] = c->docs[i][i1].paras[j].words_local[l];
				}
				if(c_tr->max_length[i] < c_tr->docs[i][i1].paras[j].total ){
					c_tr->max_length[i] = c_tr->docs[i][i1].paras[j].total;
				}
			}
		}
		if(percent != 1){
			for(i1 = 0; i1 < c_te->ndocs[i]; i1++){
				c_te->docs[i][i1].paras = (para*)malloc(sizeof(para)*c->docs[i][num_tr+i1].nparas);
				c_te->docs[i][i1].nparas = c->docs[i][num_tr+i1].nparas;
				c_te->docs[i][i1].total = c->docs[i][num_tr+i1].total; //i1 different in c_te
				c_te->totalw += c_te->docs[i][i1].total;

				for(j = 0; j < c_te->docs[i][i1].nparas; j++){
					c_te->docs[i][i1].paras[j].words = (int*)malloc(sizeof(int) * c->docs[i][num_tr+i1].paras[j].total); //i1 different in c_te
					c_te->docs[i][i1].paras[j].words_local = (int*)malloc(sizeof(int) * c->docs[i][num_tr+i1].paras[j].total); //i1 different in c_te
					c_te->docs[i][i1].paras[j].total = c->docs[i][num_tr+i1].paras[j].total;
					c_te->total[i] += c->docs[i][num_tr+i1].paras[j].total;
					for(l = 0; l < c->docs[i][num_tr+i1].paras[j].total; l++){ //i1 different in c_te
						c_te->docs[i][i1].paras[j].words[l] = c->docs[i][num_tr+i1].paras[j].words[l];
						c_te->docs[i][i1].paras[j].words_local[l] = c->docs[i][num_tr+i1].paras[j].words_local[l];
					}
					if(c_te->max_length[i] < c_te->docs[i][i1].paras[j].total ){
						c_te->max_length[i] = c_te->docs[i][i1].paras[j].total;
					}
				}
			}
		}
	}

	free_corpus(c);
}

void save_local_vocabulary(Corpus* c, vocabulary* v, char* root)
{
	FILE* fileptr;
	char fileName[BUFSIZ];
	int i, j, k;
	sprintf(fileName, "%s/local_vocabulary", root);
	mkdir(fileName, S_IRUSR | S_IWUSR | S_IXUSR);
	for(i = 0; i < c->I; ++i){
		sprintf(fileName, "%s/local_vocabulary/voc_%d.dat", root, i + 1);
		fileptr = fopen(fileName, "w");
		for(j = 0; j < c->W_i[i]; ++j){
			for(k = 0; k < v->size; ++k){
				if(v->word_map[k].id == c->l2g_vocabulary[i][j]){
					fprintf(fileptr, "%d:%s\n", j, v->word_map[k].word_str);
					break;
				}
			}
		}
		fclose(fileptr);
	}
}


/*
 * write a corpus to file
 *
 */

void write_corpus(Corpus* c, char** filename)
{
	FILE* fileptr;
	int ii, i, j, k;
	doc* d;
	para* p;
	for(ii = 0; ii < c->I; ++ii){
		fileptr = fopen(filename[ii], "w");
		for (i = 0; i < c->ndocs[ii]; i++){
			d = &(c->docs[ii][i]);
			fprintf(fileptr, "doc:%d\n", d->nparas);
			for (j = 0; j < d->nparas; j++){
				p = &(d->paras[j]);
				fprintf(fileptr, "para:%d\n", p->total);
				for (k = 0; k < p->total; k++){
					fprintf(fileptr, "word:%d\n", p->words[k]);
				}
			}
		}
		fclose(fileptr);
	}
}

/*
 * Split testing corpus to two parts
 *
 */
void split_corpus(Corpus* s, Corpus* c1, Corpus* c2)
{
	int ii, i, j, l;

	c1->I = s->I;
	c1->total = (int*)malloc(sizeof(int) * s->I);
	c1->ndocs = (int*)malloc(sizeof(int) * s->I);
	c1->max_length = (int*)malloc(sizeof(int) * s->I);
	c1->W_i = (int*)malloc(sizeof(int) * s->I);
	c1->l2g_vocabulary = (int**)malloc(sizeof(int*) * s->I);
	c1->docs = (doc**)malloc(sizeof(doc*) * s->I);
	c1->totalw = 0;
	c2->I = s->I;
	c2->total = (int*)malloc(sizeof(int) * s->I);
	c2->ndocs = (int*)malloc(sizeof(int) * s->I);
	c2->max_length = (int*)malloc(sizeof(int) * s->I);
	c2->W_i = (int*)malloc(sizeof(int) * s->I);
	c2->l2g_vocabulary = (int**)malloc(sizeof(int*) * s->I);
	c2->docs = (doc**)malloc(sizeof(doc*) * s->I);
	c2->totalw = 0;

	for(ii = 0; ii < s->I; ++ii){

		c1->docs[ii] = (doc*)malloc(sizeof(doc)*s->ndocs[ii]);
		c1->ndocs[ii] = s->ndocs[ii];
		c1->total[ii] = 0;
		c1->max_length[ii] = 0;
		c1->W_i[ii] = s->W_i[ii];
		c1->l2g_vocabulary[ii] = (int*)malloc(sizeof(int) * c1->W_i[ii]);
		for(i = 0; i < c1->W_i[ii]; ++i)
			c1->l2g_vocabulary[ii][i] = s->l2g_vocabulary[ii][i];

		c2->docs[ii] = (doc*)malloc(sizeof(doc)*s->ndocs[ii]);
		c2->ndocs[ii] = s->ndocs[ii];
		c2->total[ii] = 0;
		c2->max_length[ii] = 0;
		c2->W_i[ii] = s->W_i[ii];
		c2->l2g_vocabulary[ii] = (int*)malloc(sizeof(int) * c2->W_i[ii]);
		for(i = 0; i < c2->W_i[ii]; ++i)
			c2->l2g_vocabulary[ii][i] = s->l2g_vocabulary[ii][i];

		for(i = 0; i < s->ndocs[ii]; i++){
			c1->docs[ii][i].paras = (para*)malloc(sizeof(para)*s->docs[ii][i].nparas);
			c1->docs[ii][i].nparas = s->docs[ii][i].nparas;
			c1->docs[ii][i].total = 0;

			c2->docs[ii][i].paras = (para*)malloc(sizeof(para)*s->docs[ii][i].nparas);
			c2->docs[ii][i].nparas = s->docs[ii][i].nparas;
			c2->docs[ii][i].total = 0;

			for(j = 0; j < s->docs[ii][i].nparas; j++){
				c1->docs[ii][i].paras[j].words = (int*)malloc(sizeof(int) * 1);
				c1->docs[ii][i].paras[j].words_local = (int*)malloc(sizeof(int) * 1);
				c1->docs[ii][i].paras[j].total = 0;

				c2->docs[ii][i].paras[j].words = (int*)malloc(sizeof(int) * 1);
				c2->docs[ii][i].paras[j].words_local = (int*)malloc(sizeof(int) * 1);
				c2->docs[ii][i].paras[j].total = 0;
				for(l = 0; l < s->docs[ii][i].paras[j].total; l++){
					if(l%2 == 0){
						c1->docs[ii][i].paras[j].words = (int*)realloc(c1->docs[ii][i].paras[j].words, sizeof(int) * (c1->docs[ii][i].paras[j].total + 1));
						c1->docs[ii][i].paras[j].words_local = (int*)realloc(c1->docs[ii][i].paras[j].words_local, sizeof(int) * (c1->docs[ii][i].paras[j].total + 1));
						c1->docs[ii][i].paras[j].words[c1->docs[ii][i].paras[j].total] = s->docs[ii][i].paras[j].words[l];
						c1->docs[ii][i].paras[j].words_local[c1->docs[ii][i].paras[j].total] = s->docs[ii][i].paras[j].words_local[l];
						c1->docs[ii][i].paras[j].total++;
					}else if(l%2 == 1){
						c2->docs[ii][i].paras[j].words = (int*)realloc(c2->docs[ii][i].paras[j].words, sizeof(int) * (c2->docs[ii][i].paras[j].total + 1));
						c2->docs[ii][i].paras[j].words_local = (int*)realloc(c2->docs[ii][i].paras[j].words_local, sizeof(int) * (c2->docs[ii][i].paras[j].total + 1));
						c2->docs[ii][i].paras[j].words[c2->docs[ii][i].paras[j].total] = s->docs[ii][i].paras[j].words[l];
						c2->docs[ii][i].paras[j].words_local[c2->docs[ii][i].paras[j].total] = s->docs[ii][i].paras[j].words_local[l];
						c2->docs[ii][i].paras[j].total++;
					}
				}
				if(c1->max_length[ii] < c1->docs[ii][i].paras[j].total ){
					c1->max_length[ii] = c1->docs[ii][i].paras[j].total;
				}
				c1->docs[ii][i].total += c1->docs[ii][i].paras[j].total;

				if(c2->max_length[ii] < c2->docs[ii][i].paras[j].total ){
					c2->max_length[ii] = c2->docs[ii][i].paras[j].total ;
				}
				c2->docs[ii][i].total += c2->docs[ii][i].paras[j].total;
			}
			c1->total[ii] += c1->docs[ii][i].total;
			c1->totalw += c1->docs[ii][i].total;
			c2->total[ii] += c2->docs[ii][i].total;
			c2->totalw += c2->docs[ii][i].total;
		}
	}
}

/*
 * create word_map which map the term in vocabulary to index
 *
 * file: the file contains the corpus vocabulary in the format: index:word_string
 */
vocabulary* read_vocabulary(char* file1)
{
	char word_str[BUFSIZ];
	int word_index, count;
	FILE* fileptr;
	vocabulary* v;

	printf("Reading vocabulary ...\n");
	v = (vocabulary*) malloc(sizeof(vocabulary));
	if(file1){
		fileptr = fopen(file1, "r");
		if(!fileptr){
			printf("Cannot open vocabulary file %s\n", file1);
			exit(0);
		}
		count = 0;
		v->word_map = (id_2_word*) malloc(sizeof(id_2_word) * 1);
		while ((fscanf(fileptr, "%d:%s", &word_index, word_str)) != EOF){
			v->word_map = (id_2_word*) realloc(v->word_map, sizeof(id_2_word) * (1 + count));
			v->word_map[count].id = word_index;
			strcpy(v->word_map[count].word_str, word_str);
			assert(count == word_index);
			count++;
		}
		fclose(fileptr);
		v->size = count;
		printf("the size of the sift vocabulary: %d\n", v->size);

	}else{
		printf("vocabulary file null...\n");
		exit(0);
	}
	v_size = v->size;
	//printf("TEST FOR CRASH");

	return (v);
}

/*
 * write vocabulary to a file, be used in test
 *
 * v:  vocabulary
 * vocabulary_file: file to write v
 */
void write_vocabulary(vocabulary* v, char* file1)
{
	FILE* fileptr;
	int i;

	fileptr = fopen(file1, "w");
	for (i = 0; i < v->size; i++){
		fprintf(fileptr, "%d:%s\n", v->word_map[i].id, v->word_map[i].word_str);
	}
	fclose(fileptr);
}

void free_corpus(Corpus* c){
	int ii, i, j;
	for(ii = 0; ii < c->I; ++ii){
		for (i = 0; i < c->ndocs[ii]; i++){
			for (j = 0; j < c->docs[ii][i].nparas; j++){
				free(c->docs[ii][i].paras[j].words);
				free(c->docs[ii][i].paras[j].words_local);
			}
			free(c->docs[ii][i].paras);
		}
		free(c->docs[ii]);
		free(c->l2g_vocabulary[ii]);
	}
	free(c->W_i);
	free(c->max_length);
	free(c->ndocs);
	free(c->total);
	free(c->l2g_vocabulary);
	free(c->docs);
	free(c);
}

void free_vocabulary(vocabulary* v){
	free(v->word_map);
	free(v);
}
