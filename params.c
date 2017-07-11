/*
 * params.c
 *
 */
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include "params.h"

dpyp_params PARAMS;

void new_params(int I, int K, dpyp_params* params)
{
	int i;
	params->I = I;
	params->K_i = (int*)malloc(sizeof(int) * I);
	params->W_i = (int*)malloc(sizeof(int) * I);
	params->a = (double*)malloc(sizeof(double) * K);
	params->b = (double*)malloc(sizeof(double) * K);
	params->alpha = (double*)malloc(sizeof(double) * I);
	params->top_words = (int*)malloc(sizeof(int) * I);
	params->train_data = (char**)malloc(sizeof(char*) * I);
	for(i = 0; i < I; ++i)
		params->train_data[i] = (char*)malloc(sizeof(char) * 50);
}

void free_params(dpyp_params* params)
{
	int i;
	free(params->K_i);
	free(params->W_i);
	free(params->a);
	free(params->b);
	free(params->alpha);
	free(params->top_words);
	for(i = 0; i < params->I; ++i)
		free(params->train_data[i]);
	free(params->train_data);
}

void read_pmi_params(char* filename)
{
	FILE* fileptr;
	int tmp, i;

	fileptr = fopen(filename, "r");
	tmp = fscanf(fileptr, "#models: %d\n", &(PARAMS.I));
	assert(tmp == 1);
	tmp = fscanf(fileptr, "topk: %d\n", &(PARAMS.burn_in));
	assert(tmp == 1);
	PARAMS.train_data = (char**)malloc(sizeof(char*) * PARAMS.I);
	tmp = fscanf(fileptr, "models:");
	for(i = 0; i < PARAMS.I; ++i){
		PARAMS.train_data[i] = (char*)malloc(sizeof(char) * 100);
		tmp = fscanf(fileptr, " %s", PARAMS.train_data[i]);
		assert(tmp == 1);
	}
	fclose(fileptr);
}

void read_test_params(char* filename)
{
	FILE* fileptr;
	int tmp, i;

	fileptr = fopen(filename, "r");
	tmp = fscanf(fileptr, "I: %d\n", &(PARAMS.I));
	assert(tmp == 1);
	tmp = fscanf(fileptr, "gibbs max iter: %d\n", &(PARAMS.gibbs_max_iter));
	assert(tmp == 1);
	tmp = fscanf(fileptr, "burn in: %d\n", &(PARAMS.burn_in));
	assert(tmp == 1);
	tmp = fscanf(fileptr, "sampling lag: %d\n", &(PARAMS.sampling_lag));
	assert(tmp == 1);
	PARAMS.K_i = (int*)malloc(sizeof(int) * 1);
	PARAMS.W_i = (int*)malloc(sizeof(int) * 1);
	PARAMS.a = (double*)malloc(sizeof(double) * 1);
	PARAMS.b = (double*)malloc(sizeof(double) * 1);
	PARAMS.alpha = (double*)malloc(sizeof(double) * 1);
	PARAMS.top_words = (int*)malloc(sizeof(int) * 1);
	PARAMS.train_data = (char**)malloc(sizeof(char*) * PARAMS.I);
	tmp = fscanf(fileptr, "train files:");
	for(i = 0; i < PARAMS.I; ++i){
		PARAMS.train_data[i] = (char*)malloc(sizeof(char) * 30);
		tmp = fscanf(fileptr, " %s", PARAMS.train_data[i]);
		assert(tmp == 1);
	}
	fclose(fileptr);
}

void read_params(char* filename) {
	FILE* fileptr;
	int tmp, i;

	fileptr = fopen(filename, "r");
	tmp = fscanf(fileptr, "I: %d\n", &(PARAMS.I));
	assert(tmp == 1);
	tmp = fscanf(fileptr, "gamma: %lf\n", &(PARAMS.gamma));
	assert(tmp == 1);
	tmp = fscanf(fileptr, "train percent: %lf\n", &(PARAMS.tr_percent));
	assert(tmp == 1);
	tmp = fscanf(fileptr, "gibbs max iter: %d\n", &(PARAMS.gibbs_max_iter));
	assert(tmp == 1);
	tmp = fscanf(fileptr, "burn in: %d\n", &(PARAMS.burn_in));
	assert(tmp == 1);
	tmp = fscanf(fileptr, "sampling lag: %d\n", &(PARAMS.sampling_lag));
	assert(tmp == 1);
	PARAMS.K_i = (int*)malloc(sizeof(int) * PARAMS.I);
	tmp = fscanf(fileptr, "K_i: %d", &(PARAMS.K_i[0]));
	assert(tmp == 1);
	for(i = 1; i < PARAMS.I; ++i)
		PARAMS.K_i[i] = PARAMS.K_i[0];
	tmp = fscanf(fileptr, "\n");
	PARAMS.W_i = (int*)malloc(sizeof(int) * PARAMS.I);
	tmp = fscanf(fileptr, "W_i: %d", &(PARAMS.W_i[0]));
	assert(tmp == 1);
	for(i = 1; i < PARAMS.I; ++i)
		PARAMS.W_i[i] = PARAMS.W_i[0];
	tmp = fscanf(fileptr, "\n");
	PARAMS.a = (double*)malloc(sizeof(double) * PARAMS.K_i[0]);
	tmp = fscanf(fileptr, "a: %lf", &(PARAMS.a[0]));
	assert(tmp == 1);
	for(i = 1; i < PARAMS.K_i[0]; ++i)
		PARAMS.a[i] = PARAMS.a[0];
	tmp = fscanf(fileptr, "\n");
	PARAMS.b = (double*)malloc(sizeof(double) * PARAMS.K_i[0]);
	tmp = fscanf(fileptr, "b: %lf", &(PARAMS.b[0]));
	assert(tmp == 1);
	for(i = 1; i < PARAMS.K_i[0]; ++i)
		PARAMS.b[i] = PARAMS.b[0];
	tmp = fscanf(fileptr, "\n");
	PARAMS.alpha = (double*)malloc(sizeof(double) * PARAMS.I);
	tmp = fscanf(fileptr, "alpha: %lf", &(PARAMS.alpha[0]));
	assert(tmp == 1);
	for(i = 1; i < PARAMS.I; ++i)
		PARAMS.alpha[i] = PARAMS.alpha[0];
	tmp = fscanf(fileptr, "\n");
	PARAMS.top_words = (int*)malloc(sizeof(int) * PARAMS.I);
	tmp = fscanf(fileptr, "top words: %d", &(PARAMS.top_words[0]));
	assert(tmp == 1);
	for(i = 1; i < PARAMS.I; ++i){
		PARAMS.top_words[i] = PARAMS.top_words[0];
	}
	tmp = fscanf(fileptr, "\n");
	PARAMS.train_data = (char**)malloc(sizeof(char*) * PARAMS.I);
	tmp = fscanf(fileptr, "train files:");
	for(i = 0; i < PARAMS.I; ++i){
		PARAMS.train_data[i] = (char*)malloc(sizeof(char) * 30);
		tmp = fscanf(fileptr, " %s", PARAMS.train_data[i]);
		assert(tmp == 1);
	}
	fclose(fileptr);
}

void write_params(char* filename) {
	FILE* fileptr;
	int i;
	fileptr = fopen(filename, "w");

	fprintf(fileptr, "I: %d\n", PARAMS.I);
	fprintf(fileptr, "gamma: %lf\n", PARAMS.gamma);
	fprintf(fileptr, "train percent: %lf\n", PARAMS.tr_percent);
	fprintf(fileptr, "gibbs max iter: %d\n", PARAMS.gibbs_max_iter);
	fprintf(fileptr, "burn in: %d\n", PARAMS.burn_in);
	fprintf(fileptr, "sampling lag: %d\n", PARAMS.sampling_lag);
	fprintf(fileptr, "K_i:");
	for(i = 0; i < PARAMS.I; ++i)
		fprintf(fileptr, " %d", PARAMS.K_i[i]);
	fprintf(fileptr, "\n");
	fprintf(fileptr, "W_i:");
	for(i = 0; i < PARAMS.I; ++i)
		fprintf(fileptr, " %d", PARAMS.W_i[i]);
	fprintf(fileptr, "\n");
	fprintf(fileptr, "a:");
	for(i = 0; i < PARAMS.K_i[0]; ++i)
		fprintf(fileptr, " %lf", PARAMS.a[i]);
	fprintf(fileptr, "\n");
	fprintf(fileptr, "b:");
	for(i = 0; i < PARAMS.K_i[0]; ++i)
		fprintf(fileptr, " %lf", PARAMS.b[i]);
	fprintf(fileptr, "\n");
	fprintf(fileptr, "alpha:");
	for(i = 0; i < PARAMS.I; ++i)
		fprintf(fileptr, " %lf", PARAMS.alpha[i]);
	fprintf(fileptr, "\n");
	fprintf(fileptr, "top words:");
	for(i = 0; i < PARAMS.I; ++i)
		fprintf(fileptr, " %d", PARAMS.top_words[i]);
	fprintf(fileptr, "\n");
	fprintf(fileptr, "train files:");
	for(i = 0; i < PARAMS.I; ++i){
		fprintf(fileptr, " %s", PARAMS.train_data[i]);
	}
	fprintf(fileptr, "\n");
	fclose(fileptr);
}

void print_params() {
	int i;
	printf(">>>>>> parameters <<<<<<\n");

	printf("I: %d\n", PARAMS.I);
	printf("gamma: %lf\n", PARAMS.gamma);
	printf("train percent: %lf\n", PARAMS.tr_percent);
	printf("gibbs max iter: %d\n", PARAMS.gibbs_max_iter);
	printf("burn in: %d\n", PARAMS.burn_in);
	printf("sampling lag: %d\n", PARAMS.sampling_lag);
	printf("K_i:");
	for(i = 0; i < PARAMS.I; ++i)
		printf(" %d", PARAMS.K_i[i]);
	printf("\n");
	printf("W_i:");
	for(i = 0; i < PARAMS.I; ++i)
		printf(" %d", PARAMS.W_i[i]);
	printf("\n");
	printf("a:");
	for(i = 0; i < PARAMS.K_i[0]; ++i)
		printf(" %lf", PARAMS.a[i]);
	printf("\n");
	printf("b:");
	for(i = 0; i < PARAMS.K_i[0]; ++i)
		printf(" %lf", PARAMS.b[i]);
	printf("\n");
	printf("alpha:");
	for(i = 0; i < PARAMS.I; ++i)
		printf(" %lf", PARAMS.alpha[i]);
	printf("\n");
	printf("top words:");
	for(i = 0; i < PARAMS.I; ++i)
		printf(" %d", PARAMS.top_words[i]);
	printf("\n");
	printf("train files:");
	for(i = 0; i < PARAMS.I; ++i){
		printf(" %s", PARAMS.train_data[i]);
	}
	printf("\n");

	printf(">>>>>> parameters <<<<<<\n");
}

void set_paras(int I, int *K_i, int* W_i, double *a, double *b, double *alpha, double gamma,
		double tr_per, int gibbs_max_iter, int burn_in, int sampling_lag, int *top_words)
{
	int i;
	PARAMS.I = I;
	PARAMS.gamma = gamma;
	PARAMS.tr_percent = tr_per;
	PARAMS.gibbs_max_iter = gibbs_max_iter;
	PARAMS.burn_in = burn_in;
	PARAMS.sampling_lag = sampling_lag;
	for(i = 0; i < K_i[0]; ++i){
		PARAMS.a[i] = a[i];
		PARAMS.b[i] = b[i];
	}
	for(i = 0; i < I; ++i){
		PARAMS.K_i[i] = K_i[i];
		PARAMS.W_i[i] = W_i[i];
		PARAMS.alpha[i] = alpha[i];
		PARAMS.top_words[i] = top_words[i];
	}
}

void default_params() {

}
