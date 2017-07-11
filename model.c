/*
 * model.c
 *
 */
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <gsl/gsl_sf.h>
#include <gsl/gsl_math.h>
#include <unistd.h>

#include "model.h"
#include "util.h"
#include "hash.h"

Assignment* ass_tmp;

/*
 * k:	number of topics
 * v: 	size of vocabulary
 */
Model* new_model(int I, int v, int *K_i, int* W_i, int inf)
{
	int i;
	Model* model = (Model*) malloc(sizeof(Model));
	model->I = I;
	model->v = v;
	model->K_I = (int*)malloc(sizeof(int) * I);
	model->W_I = (int*)malloc(sizeof(int) * I);
	model->a = (double*)malloc(sizeof(double) * K_i[0]);
	model->b = (double*)malloc(sizeof(double) * K_i[0]);
	if(inf == 1){
		model->rest_P = NULL;
	}else{
		model->rest_P = gsl_matrix_calloc(I, v);
	}
	for(i = 0; i < I; ++i){
		model->K_I[i] = K_i[i];
		model->W_I[i] = W_i[i];
	}
	model->alpha = (gsl_vector**)malloc(sizeof(gsl_vector*) * I);
	for(i = 0; i < I; ++i)
		model->alpha[i] = gsl_vector_calloc(K_i[i]);
	model->alpha_sum = (double*)malloc(sizeof(double) * I);
	model->gamma = gsl_vector_calloc(v);
	model->phi = gsl_matrix_calloc(K_i[0], v);
	model->phi_i = (gsl_matrix**)malloc(sizeof(gsl_matrix*) * I);
	if(inf != 1){
		model->P_k = (gsl_vector**)malloc(sizeof(gsl_vector*) * I);
	}else{
		model->P_k = NULL;
	}
	for(i = 0; i < I; ++i){
		if(inf == 1)
			model->phi_i[i] = gsl_matrix_calloc(K_i[i], v);
		else
			model->phi_i[i] = gsl_matrix_calloc(K_i[i], W_i[i]);
		if(inf != 1){
			model->P_k[i] = gsl_vector_calloc(W_i[i]);
		}
	}
	if(inf == 1){
		model->P = NULL;
		model->P_ind = NULL;
		model->P_size = NULL;
		model->P_size_all = NULL;
		model->q = NULL;
		model->vL = NULL;
		model->P_all = NULL;
		model->P_idx_all = NULL;
		model->P_idx_all = NULL;
	}else{
		model->P = (double***)malloc(sizeof(double**) * I);
		model->P_ind = (int***)malloc(sizeof(int**) * I);
		model->P_size = (int**)malloc(sizeof(int*) * I);
		for(i = 0; i < I; ++i){
			model->P[i] = (double**)malloc(sizeof(double*) * W_i[i]);
			model->P_ind[i] = (int**)malloc(sizeof(int*) * W_i[i]);
			model->P_size[i] = (int*)malloc(sizeof(int) * W_i[i]);
		}
		int k;
		model->q = (int****)malloc(sizeof(int***) * I);
		model->vL = (vlist****)malloc(sizeof(vlist***) * I);
		for(i = 0; i < I; ++i){
			model->q[i] = (int***)malloc(sizeof(int**) * K_i[i]);
			model->vL[i] = (vlist***)malloc(sizeof(vlist**) * K_i[i]);
			for(k = 0; k < K_i[i]; ++k){
				model->q[i][k] = (int**)malloc(sizeof(int*) * W_i[i]);
				model->vL[i][k] = (vlist**)malloc(sizeof(vlist*) * W_i[i]);
			}
		}
		model->P_all = (double**)malloc(sizeof(double*) * model->v);
		model->P_idx_all = (int**)malloc(sizeof(int*) * model->v);
		model->P_size_all = (int*)malloc(sizeof(int) * model->v);
		for(i = 0; i < model->v; ++i){
			model->P_all[i] = NULL;
			model->P_idx_all[i] = NULL;
		}
	}

	return model;
}
/*
 * Initialize all count tables for statistics
 *

 */
Cts* new_cts(int I, int * K_i, int *W_i, int v, Corpus* c)
{
	int i, j;

	Cts* cts = (Cts*) malloc(sizeof(Cts));
	cts->B = (double*)malloc(sizeof(double) * K_i[0]);
	cts->sum_T = (double*)malloc(sizeof(double) * K_i[0]);
	cts->T_KQ = (double**)malloc(sizeof(double*) * K_i[0]);
	for(i = 0; i < K_i[0]; ++i){
		cts->T_KQ[i] = (double*)malloc(sizeof(double) * v);
	}
	cts->n = (int***)malloc(sizeof(int**) * I);
	cts->N = (int**)malloc(sizeof(int*) * I);
	for(i = 0; i < I; ++i){
		cts->n[i] = (int**)malloc(sizeof(int*) * c->ndocs[i]);
		cts->N[i] = (int*)malloc(sizeof(int) * c->ndocs[i]);
		for(j = 0; j < c->ndocs[i]; ++j){
			cts->n[i][j] = (int*)malloc(sizeof(int) * K_i[i]);
		}
	}
	cts->top_cts = (Top_cts*)malloc(sizeof(Top_cts) * I);
	for(i = 0; i < I; ++i){
		cts->top_cts[i].m = (int**)malloc(sizeof(int*) * K_i[i]);
		cts->top_cts[i].M = (int*)malloc(sizeof(int) * K_i[i]);
		cts->top_cts[i].t = (int**)malloc(sizeof(int*) * K_i[i]);
		cts->top_cts[i].T_K = (int*)malloc(sizeof(int) * K_i[i]);
		cts->top_cts[i].T_W = (int*)malloc(sizeof(int) * W_i[i]);
		for(j = 0; j < K_i[i]; ++j){
			cts->top_cts[i].m[j] = (int*)malloc(sizeof(int) * W_i[i]);
			cts->top_cts[i].t[j] = (int*)malloc(sizeof(int) * W_i[i]);
		}
	}
	return cts;
}

/*
 * Initialize the topic assignment of each word-token
 *
 */
Assignment* new_assignment(Corpus* c)
{
	int i, j;

	Assignment* ass = (Assignment*) malloc(sizeof(Assignment));
	ass->topic_ass = (int***) malloc(sizeof(int**) * c->I);
	for (i = 0; i < c->I; i++){
		ass->topic_ass[i] = (int**) malloc(sizeof(int*) * c->ndocs[i]);
		for (j = 0; j < c->ndocs[i]; j++){
			ass->topic_ass[i][j] = (int*)malloc(sizeof(int) * c->docs[i][j].paras[0].total);
		}
	}
	ass_tmp = ass;
	return ass;
}

/*
 * Initialzie the estimator for compute mu, nu, and phi
 *
 */
Estimator* new_estimator(Corpus* c, int I, int* K_i, int *W_i, int v, int do_mu, int do_phi)
{
	int i;

	Estimator* est = (Estimator*)malloc(sizeof(Estimator));

	if(do_mu){
		est->mu = (gsl_matrix**)malloc(sizeof(gsl_matrix*) * I);
		for(i = 0; i < I; ++i)
			est->mu[i] = gsl_matrix_calloc(c->ndocs[i], K_i[i]);
	}else{
		est->mu = NULL;
	}

	if(do_phi){
		est->phi_i = (gsl_matrix**)malloc(sizeof(gsl_matrix*) * I);
		for(i = 0; i < I; ++i)
			est->phi_i[i] = gsl_matrix_calloc(K_i[i], W_i[i]);
		est->phi = gsl_matrix_calloc(K_i[0], v); //gsl_matrix_calloc(tot_t, v);
	}else{
		est->phi = NULL;
		est->phi_i = NULL;
	}

	return est;
}
/*
 * randomly initialize the model
 *
 */

Model* random_init(int I, int *K_i, int *W_i, int v, double *alpha, double gamma,
		double *a, double *b, int** l2g_vocabulary, char* dir_p)
{
	int i, ii;
	double ss;
	char file[BUFSIZ];
	Model* model;
	model = new_model(I, v, K_i, W_i, 0);
	for(i = 0; i < K_i[0]; ++i){
		model->a[i] = a[i];
		model->b[i] = b[i];
	}
	gsl_matrix_set_all(model->rest_P, 1);
	if(dir_p){
		sprintf(file, "%s/P_size.txt", dir_p);
		load_P_size(file, model);
		sprintf(file, "%s/P_size_all.txt", dir_p);
		load_P_size_all(file, model);
	}
	for(i = 0; i < I; ++i){
		model->alpha_sum[i] = alpha[i] * K_i[i];
		model->gamma_sum = gamma * v;
		gsl_vector_set_all(model->alpha[i], alpha[i]);
		gsl_vector_set_all(model->gamma, gamma);

		if(dir_p){
			int i1;

			for(i1 = 0; i1 < W_i[i]; ++i1){
				model->P[i][i1] = (double*)malloc(sizeof(double) * model->P_size[i][i1]);
				model->P_ind[i][i1] = (int*)malloc(sizeof(int) * model->P_size[i][i1]);
			}
			sprintf(file, "%s/P_ind_%d.txt", dir_p, i);
			load_P_ind(file, model, i);
			sprintf(file, "%s/P_K_%d.txt", dir_p, i);
			load_P(file, model, i);

		}else{
			int i1;
			for(i1 = 0; i1 < W_i[i]; ++i1){
				model->P[i][i1] = (double*)malloc(sizeof(double) * 1);
				model->P_ind[i][i1] = (int*)malloc(sizeof(int) * 1);
				model->P_size[i][i1] = 1;
				model->P[i][i1][0] =  1;
				model->P_ind[i][i1][0] = l2g_vocabulary[i][i1];
			}
		}
		int k, w, vv;
		for(k = 0; k < K_i[i]; ++k){
			for(w = 0; w < W_i[i]; ++w){
				model->vL[i][k][w] = NULL;
				model->q[i][k][w] = (int*)malloc(sizeof(int) * model->P_size[i][w]);
				for(vv = 0; vv < model->P_size[i][w]; ++vv)
					model->q[i][k][w][vv] = 0;
			}
		}
		for(ii = 0; ii < W_i[i]; ++ii){
			int j;
			ss = 0;
			for(j = 0; j < model->P_size[i][ii]; ++j){
				ss += model->P[i][ii][j];
				minc(model->rest_P, i, model->P_ind[i][ii][j], -model->P[i][ii][j]);
			}
			vset(model->P_k[i], ii, ss);
		}
	}

	if(dir_p){
		int i1;
		for(i1 = 0; i1 < model->v; ++i1){
			model->P_all[i1] = (double*)malloc(sizeof(double) * model->P_size_all[i1]);
			model->P_idx_all[i1] = (int*)malloc(sizeof(int) * model->P_size_all[i1]);
		}

		sprintf(file, "%s/P_idx_all.txt", dir_p);
		load_P_ind_all(file, model);
		sprintf(file, "%s/P_K_all.txt", dir_p);
		load_P_all(file, model);
	}else{
		int i1;
		for(i1 = 0; i1 < model->v; ++i1){
			model->P_all[i1] = (double*)malloc(sizeof(double) * 1);
			model->P_idx_all[i1] = (int*)malloc(sizeof(int) * 1);
			model->P_size_all[i1] = 1;
			model->P_all[i1][0] =  1;
			model->P_idx_all[i1][0] = i1;
		}
	}
	return model;
}

void save_P(char* file, Model* model, int i)
{
	int w, v;
	FILE* fileptr = fopen(file, "w");
	for(w = 0; w < model->W_I[i]; ++w){
		for(v = 0; v < model->P_size[i][w]; ++v){
			fprintf(fileptr, "%lf ", model->P[i][w][v]);
		}
		fprintf(fileptr, "\n");
	}
	fclose(fileptr);
}

void load_P(char* file, Model* model, int i)
{
	int w, v, tmp;
	FILE* fileptr = fopen(file, "r");
	for(w = 0; w < model->W_I[i]; ++w){
		for(v = 0; v < model->P_size[i][w]; ++v){
			tmp = fscanf(fileptr, "%lf ", &(model->P[i][w][v]));
			assert(tmp == 1);
		}
		tmp = fscanf(fileptr, "\n");
	}
	fclose(fileptr);
}

void load_P_all(char* file, Model* model)
{
	int w, v, tmp;
	FILE* fileptr = fopen(file, "r");
	for(w = 0; w < model->v; ++w){
		for(v = 0; v < model->P_size_all[w]; ++v){
			tmp = fscanf(fileptr, "%lf ", &(model->P_all[w][v]));
			assert(tmp == 1);
		}
		tmp = fscanf(fileptr, "\n");
	}
	fclose(fileptr);
}

void save_P_ind(char* file, Model* model, int i)
{
	int w, v;
	FILE* fileptr = fopen(file, "w");
	for(w = 0; w < model->W_I[i]; ++w){
		for(v = 0; v < model->P_size[i][w]; ++v){
			fprintf(fileptr, "%d ", model->P_ind[i][w][v]);
		}
		fprintf(fileptr, "\n");
	}
	fclose(fileptr);
}

void load_P_ind(char* file, Model* model, int i)
{
	int w, v, tmp;
	FILE* fileptr = fopen(file, "r");
	for(w = 0; w < model->W_I[i]; ++w){
		for(v = 0; v < model->P_size[i][w]; ++v){
			tmp = fscanf(fileptr, "%d ", &(model->P_ind[i][w][v]));
			assert(tmp == 1);
		}
		tmp = fscanf(fileptr, "\n");
	}
	fclose(fileptr);
}

void load_P_ind_all(char* file, Model* model)
{
	int w, v, tmp;
	FILE* fileptr = fopen(file, "r");
	for(w = 0; w < model->v; ++w){
		for(v = 0; v < model->P_size_all[w]; ++v){
			tmp = fscanf(fileptr, "%d ", &(model->P_idx_all[w][v]));
			assert(tmp == 1);
		}
		tmp = fscanf(fileptr, "\n");
	}
	fclose(fileptr);
}

void save_P_size(char* file, Model* model)
{
	int i, w;
	FILE* fileptr = fopen(file, "w");
	for(i = 0; i < model->I; ++i){
		for(w = 0; w < model->W_I[i]; ++w){
			fprintf(fileptr, "%d ", model->P_size[i][w]);
		}
		fprintf(fileptr, "\n");
	}
	fclose(fileptr);
}

void load_P_size(char* file, Model* model)
{
	int i, w, tmp;
	FILE* fileptr = fopen(file, "r");
	for(i = 0; i < model->I; ++i){
		for(w = 0; w < model->W_I[i]; ++w){
			tmp = fscanf(fileptr, "%d ", &(model->P_size[i][w]));
			assert(tmp == 1);
		}
		tmp = fscanf(fileptr, "\n");
	}
	fclose(fileptr);
}

void load_P_size_all(char* file, Model* model)
{
	int v, tmp;
	FILE* fileptr = fopen(file, "r");
	for(v = 0; v < model->v; ++v){
		tmp = fscanf(fileptr, "%d ", &(model->P_size_all[v]));
		assert(tmp == 1);
	}
	tmp = fscanf(fileptr, "\n");
	fclose(fileptr);
}

/*
 * read model
 *
 */
Model* read_model(char* file)
{
	char str[BUFSIZ];
	FILE* fileptr;
	Model* model;
	int v, i, I, *K_i, *W_i;
	int tmp;
	double alp, gamma_sum;

	sprintf(str, "%s/model_other.txt", file);
	fileptr = fopen(str, "r");
	if (!fileptr){
		printf("File %s/model_other.txt does not exist\n", file);
		exit(0);
	}
	tmp = fscanf(fileptr, "model->I: %d\n", &I);
	assert(tmp == 1);

	K_i = (int*)malloc(sizeof(int) * I);
	W_i = (int*)malloc(sizeof(int) * I);

	tmp = fscanf(fileptr, "model->v: %d\n", &v);
	assert(tmp == 1);
	tmp = fscanf(fileptr, "model->gamma_sum: %lf\n", &gamma_sum);
	assert(tmp == 1);
	tmp = fscanf(fileptr, "model->K_i:");
	for(i = 0; i < I; ++i){
		tmp = fscanf(fileptr, " %d", &K_i[i]);
		assert(tmp == 1);
	}
	tmp = fscanf(fileptr, "\n");
	tmp = fscanf(fileptr, "model->W_i:");
	for(i = 0; i < I; ++i){
		tmp = fscanf(fileptr, " %d", &W_i[i]);
		assert(tmp == 1);
	}
	tmp = fscanf(fileptr, "\n");
	model = new_model(I, v, K_i, W_i, 1);

	model->rest_P = NULL;

	model->gamma_sum = gamma_sum;

	tmp = fscanf(fileptr, "model->a:");
	for(i = 0; i < K_i[0]; ++i){
		tmp = fscanf(fileptr, " %lf", &model->a[i]);
		assert(tmp == 1);
	}
	tmp = fscanf(fileptr, "\n");
	tmp = fscanf(fileptr, "model->b:");
	for(i = 0; i < K_i[0]; ++i){
		tmp = fscanf(fileptr, " %lf", &model->b[i]);
		assert(tmp == 1);
	}
	tmp = fscanf(fileptr, "\n");
	tmp = fscanf(fileptr, "model->alpha_sum:");
	for(i = 0; i < I; ++i){
		tmp = fscanf(fileptr, " %lf", &alp);
		assert(tmp == 1);
		model->alpha_sum[i] = alp;
	}
	tmp = fscanf(fileptr, "\n");
	fclose(fileptr);

	sprintf(str, "%s/model_gamma.txt", file);
	scanf_vector(str, model->gamma);
	sprintf(str, "%s/model_phi.txt", file);
	scanf_matrix(str, model->phi);
	//sprintf(str, "%s/model_P_size.txt", file);
	//load_P_size(str, model);
	for(i = 0; i < I; ++i){
		sprintf(str, "%s/model_alpha_%d.txt", file, i);
		scanf_vector(str, model->alpha[i]);

		//sprintf(str, "%s/model_P_k_%d.txt", file, i);
		//scanf_vector(str, model->P_k[i]);

		sprintf(str, "%s/model_phi_%d.txt", file, i);
		scanf_matrix(str, model->phi_i[i]);

		int w;
		int ii;
		double ss;
		for(ii = 0; ii < model->K_I[0]; ++ii){
			ss = 0;
			for(w = 0; w < model->v; ++w){
				ss += mget(model->phi_i[i], ii, w);
			}
			for(w = 0; w < model->v; ++w){
				mset(model->phi_i[i], ii, w, mget(model->phi_i[i], ii, w)/ss);
			}
		}

		/*for(w = 0; w < model->W_I[i]; ++w){
			model->P[i][w] = (double*)malloc(sizeof(double) * model->P_size[i][w]);
			model->P_ind[i][w] = (int*)malloc(sizeof(int) * model->P_size[i][w]);
		}
		sprintf(str, "%s/model_P_%d.txt", file, i);
		load_P(str, model, i);
		sprintf(str, "%s/model_P_ind_%d.txt", file, i);
		load_P_ind(str, model, i);*/
	}

	return model;
}
/*
 * function: to save model
 */
void save_model(Model* model, Corpus* c, Cts* cts, gsl_matrix* mphi, gsl_matrix** mphi_i, char* file)
{
	char str[BUFSIZ];
	FILE* fileptr;
	int i;
	gsl_matrix* phi_tmp = gsl_matrix_calloc(mphi->size1, mphi->size2);

	sprintf(str, "%s/model_other.txt", file);
	fileptr = fopen(str, "w");
	fprintf(fileptr, "model->I: %d\n", model->I);
	fprintf(fileptr, "model->v: %d\n", model->v);
	fprintf(fileptr, "model->gamma_sum: %lf\n", model->gamma_sum);
	fprintf(fileptr, "model->K_i:");
	for(i = 0; i < model->I; ++i)
		fprintf(fileptr, " %d", model->K_I[i]);
	fprintf(fileptr, "\n");
	fprintf(fileptr, "model->W_i:");
	for(i = 0; i < model->I; ++i)
		fprintf(fileptr, " %d", model->W_I[i]);
	fprintf(fileptr, "\n");
	fprintf(fileptr, "model->a:");
	for(i = 0; i < model->K_I[0]; ++i)
		fprintf(fileptr, " %lf", model->a[i]);
	fprintf(fileptr, "\n");
	fprintf(fileptr, "model->b:");
	for(i = 0; i < model->K_I[0]; ++i)
		fprintf(fileptr, " %lf", model->b[i]);
	fprintf(fileptr, "\n");
	fprintf(fileptr, "model->alpha_sum:");
	for(i = 0; i < model->I; ++i)
		fprintf(fileptr, " %lf", model->alpha_sum[i]);
	fprintf(fileptr, "\n");
	fclose(fileptr);

	sprintf(str, "%s/model_gamma.txt", file);
	save_vector(str, model->gamma);
	sprintf(str, "%s/model_phi.txt", file);
	save_matrix(str, mphi);

	sprintf(str, "%s/model_P_size.txt", file);
	save_P_size(str, model);
	for(i = 0; i < model->I; ++i){
		int k1, w1;
		sprintf(str, "%s/model_phi_%d.txt", file, i);
		gsl_matrix_set_all(phi_tmp, 0);
		for(k1 = 0; k1 < phi_tmp->size1; ++k1){
			if(model->v != model->W_I[i]){
				for(w1 = 0; w1 < phi_tmp->size2; ++w1){
					double ss = 0;
					int i1;
					for(i1 = 0; i1 < model->P_size_all[w1]; ++i1){
						ss += mget(mphi, k1, model->P_idx_all[w1][i1]) * model->P_all[w1][i1];
					}
					gsl_matrix_set(phi_tmp, k1, w1,
							(model->b[k1] + cts->top_cts[i].T_K[k1] * model->a[k1]) * ss
							/ (model->b[k1] + cts->top_cts[i].M[k1]));
				}
			}
			for(w1 = 0; w1 < mphi_i[i]->size2; ++w1){
				//mset(phi_tmp, k1, c->l2g_vocabulary[i][w1], mget(mphi_i[i], k1, w1));
				minc(phi_tmp, k1, c->l2g_vocabulary[i][w1], mget(mphi_i[i], k1, w1));
			}
			double summ = 0;
			for(w1 = 0; w1 < phi_tmp->size2; ++w1){
				summ += mget(phi_tmp, k1, w1);
			}
			for(w1 = 0; w1 < phi_tmp->size2; ++w1){
				mset(phi_tmp, k1, w1, mget(phi_tmp, k1, w1) / summ);
			}
		}
		save_matrix(str, phi_tmp);

		sprintf(str, "%s/model_phi_sub%d.txt", file, i);
		save_matrix(str, mphi_i[i]);

		sprintf(str, "%s/model_alpha_%d.txt", file, i);
		save_vector(str, model->alpha[i]);

		sprintf(str, "%s/model_P_k_%d.txt", file, i);
		save_vector(str, model->P_k[i]);

		sprintf(str, "%s/model_P_%d.txt", file, i);
		save_P(str, model, i);
		sprintf(str, "%s/model_P_ind_%d.txt", file, i);
		save_P_ind(str, model, i);
	}
	gsl_matrix_free(phi_tmp);
}

void save_customer_table(Cts* cts, int I, int K, int v, char* file)
{
	int i, k, w;
	FILE *fileptr;
	char str[BUFSIZ];
	for(i = 0; i < I; i++){
		sprintf(str, "%s/model_customer_%d.txt", file, i);
		fileptr = fopen(str, "w");
		for(k = 0; k < K; k++){
			for(w = 0; w < v; w++){
				fprintf(fileptr, "%d ", cts->top_cts[i].m[k][w]);
			}
			fprintf(fileptr, "\n");
		}
		fclose(fileptr);

		sprintf(str, "%s/model_table_%d.txt", file, i);
		fileptr = fopen(str, "w");
		for(k = 0; k < K; k++){
			for(w = 0; w < v; w++){
				fprintf(fileptr, "%d ", cts->top_cts[i].t[k][w]);
			}
			fprintf(fileptr, "\n");
		}
		fclose(fileptr);
	}
}

/*
 * read and save estimator
 */
void save_estimator(Estimator* est, Corpus* c, char* file, int do_mu, int do_phi)
{
	int i;
	char str[BUFSIZ];

	if(do_mu){
		for(i = 0; i < c->I; ++i){
			sprintf(str, "%s/mu_%d.txt", file, i);
			save_matrix(str, est->mu[i]);
		}
	}

	if(do_phi){
		sprintf(str, "%s/phi.txt", file);
		save_matrix(str, est->phi);
		for(i = 0; i < c->I; ++i){
			sprintf(str, "%s/phi_%d.txt", file, i);
			save_matrix(str, est->phi_i[i]);
		}
	}
}

void printf_estimator(Estimator* est, Corpus* c, char* file, int do_mu, int do_phi)
{
	int i;
	char str[BUFSIZ];

	if(do_mu){
		for(i = 0; i < c->I; ++i){
			sprintf(str, "%s/mu_%d.txt", file, i);
			printf_matrix(str, est->mu[i]);
		}
	}

	if(do_phi){
		sprintf(str, "%s/phi.txt", file);
		printf_matrix(str, est->phi);
		for(i = 0; i < c->I; ++i){
			sprintf(str, "%s/phi_%d.txt", file, i);
			save_matrix(str, est->phi_i[i]);
		}
	}
}

Estimator* read_estimator(Corpus* c, int v, int *K_i, int *W_i, char* file, int do_mu, int do_phi)
{
	int i;
	char str[BUFSIZ];

	Estimator* est =  new_estimator(c, c->I, K_i, W_i, v, do_mu, do_phi);

	if(do_mu){
		for(i = 0; i < c->I; ++i){
			sprintf(str, "%s/mu_%d.txt", file, i);
			scanf_matrix(str, est->mu[i]);
		}
	}

	if(do_phi){
		sprintf(str, "%s/phi.txt", file);
		scanf_matrix(str, est->phi);
		for(i = 0; i < c->I; ++i){
			sprintf(str, "%s/phi_%d.txt", file, i);
			scanf_matrix(str, est->phi_i[i]);
		}
	}
	return est;
}

/*
 * read topic assignment
 */
Assignment* read_topic_assignmnet(Corpus* c, char* file)
{
	Assignment* ass = new_assignment(c);
	FILE* fileptr;
	int i, j, l, tmp;
	char str[BUFSIZ];

	printf("Reading topic assignment ......\n");
	for(i = 0; i < c->I; i++){
		sprintf(str, "%s/z_%d.txt", file, i);
		fileptr = fopen(str, "r");
		if(!fileptr){
			printf("Cannot open file %s\n", str);
			exit(0);
		}
		for(j = 0; j < c->ndocs[i]; j++){
			for(l = 0; l < c->docs[i][j].total; ++l){
				tmp = fscanf(fileptr, "%d ", &ass->topic_ass[i][j][l]);
				assert(tmp == 1);
			}
			tmp = fscanf(fileptr, "\n");
		}
		fclose(fileptr);
	}
	printf("Finished!!!\n");
	return ass;
}

void save_topic_assignmnet(Corpus* c, Assignment* ass, char* file)
{
	FILE* fileptr;
	int i, j, l;
	char str[BUFSIZ];

	printf("Saving topic assignment ......\n");
	for(i = 0; i < c->I; i++){
		sprintf(str, "%s/z_%d.txt", file, i);
		fileptr = fopen(str, "w");
		for(j = 0; j < c->ndocs[i]; j++){
			for(l = 0; l < c->docs[i][j].total; ++l){
				fprintf(fileptr, "%d ", ass->topic_ass[i][j][l]);
			}
			fprintf(fileptr, "\n");
		}
		fclose(fileptr);
	}
	printf("Finished!!!\n");
}

void print_word_asso(Model* model, vocabulary* voc, char* filename)
{
	int i, k, w, v, vv;
	FILE* fid = fopen(filename, "w");
	for(k = 0; k < model->K_I[0]; k++){
		fprintf(fid, "word association for topic %d\n", k + 1);
		for(v = 0; v < voc->size; v++){
			fprintf(fid, "--------------%s------------\n", voc->word_map[v].word_str);
			for(i = 0; i < model->I; i++){
				fprintf(fid, "group %d: ", i + 1);
				for(w = 0; w < model->W_I[i]; w++){
					for(vv = 0; vv < model->P_size[i][w]; vv++){
						if(model->P_ind[i][w][vv] == v){
							break;
						}
					}
					if(vv < model->P_size[i][w] && model->q[i][k][w][vv] > 0){
						fprintf(fid, "%s(%d) ", voc->word_map[w].word_str, model->q[i][k][w][vv]);
					}
				}
				fprintf(fid, "\n");
			}
			fprintf(fid, "\n");
		}
		fprintf(fid, "\n");
	}
	fclose(fid);
}


void print_top_words(int num_words, int begin_k, int end_k, gsl_matrix* M, Cts* cts, vocabulary* v, char* file)
{
	FILE* fileptr;
	int k, i, j;
	int word_id;
	gsl_vector* vector;
	int V, *counts, total_c = 0;
	V = v->size;

	counts = (int*)malloc(sizeof(int) * (end_k - begin_k + 1));
	for(k = 0; k < end_k - begin_k + 1; k++){
		counts[k] = 0;
	}
	for(i = 0; i < cts->I; i++){
		for(k = begin_k; k <= end_k; k++){
			counts[k - begin_k] += cts->top_cts[i].M[k];
			total_c += cts->top_cts[i].M[k];
		}
	}

	fileptr = fopen(file, "w");
	vector = gsl_vector_alloc(V);
	for (k = begin_k; k <= end_k; k++){
		fprintf(fileptr, "------topic_%d (%lf)------\n", k - begin_k, ((double)counts[k - begin_k]) / total_c);
		gsl_matrix_get_row(vector, M, k);
		sort(vector);
		for (i = 0; i < num_words; i++){
			word_id = (int) vget(vector, i);
			for (j = 0; j < V; j++){
				if (v->word_map[j].id == word_id) {
					fprintf(fileptr, "%s(%lf) ", v->word_map[j].word_str, mget(M, k, word_id));
					break;
				}
			}
		}
		fprintf(fileptr, "\n");
	}
	fclose(fileptr);
	gsl_vector_free(vector);
	free(counts);
}

void print_local_top_words(int* num_words, int begin_k, int end_k, gsl_matrix**M_i,
		Cts* cts, Corpus* c, vocabulary* v, char* file)
{
	FILE* fileptr;
	int k, i, j, d;
	int word_id;
	gsl_vector* vector;
	char str[BUFSIZ];
	int V, total_c;
	V = v->size;

	for(d = 0; d < c->I; ++d){
		total_c = 0;
		for(k = begin_k; k <= end_k; k++){
			total_c += cts->top_cts[d].M[k];
		}
		sprintf(str, "%s_sub%d.txt", file, d);
		fileptr = fopen(str, "w");
		vector = gsl_vector_alloc(M_i[d]->size2);
		for (k = begin_k; k <= end_k; k++){
			fprintf(fileptr, "------topic_%d(%lf)------\n", k - begin_k, ((double)cts->top_cts[d].M[k]) / total_c);
			gsl_matrix_get_row(vector, M_i[d], k);
			sort(vector);
			for (i = 0; i < num_words[d]; i++){
				word_id = (int) vget(vector, i);
				for (j = 0; j < V; j++){
					if (v->word_map[j].id == c->l2g_vocabulary[d][word_id]) {
						fprintf(fileptr, "%s(%lf) ", v->word_map[j].word_str, mget(M_i[d], k, word_id));
						break;
					}
				}
			}
			fprintf(fileptr, "\n");
		}
		fclose(fileptr);
		gsl_vector_free(vector);
	}
}

void calPMI(char* dir, char** model_dir, int numM, int topk, char* pmifile)
{
	int d, i, j, k, totV;
	double value;
	char s_tmp[BUFSIZ];

	FILE *fr2 = fopen(pmifile,"r");
	if (!fr2){
		printf("pmifile '%s' not read\n", pmifile);
		exit(0);
	}
	int count = 0;
	char * line = NULL;
	size_t len = 0;
	ssize_t read;
	printf("Reading PMI file %s...\n", pmifile);
	while ((read = getline(&line, &len, fr2)) != EOF) {
		 count++;
	}
	printf("Finish reading PMI file, lines = %d\n", count);
	if (line)
		 free(line);
	fclose(fr2);
	fr2 = fopen(pmifile,"r");
	if (!fr2){
		printf("pmifile '%s' not read\n", pmifile);
		exit(0);
	}
	hash_table_init();
	printf("Reading PMI file...\n");
	while (fscanf(fr2, "%d", &i) != EOF ) {
		if (fscanf(fr2, "%d", &j) != EOF) {
			if (fscanf(fr2,"%lf",&value) !=EOF) {
				sprintf(s_tmp, "%d,%d", i-1, j-1);
				hash_table_insert(s_tmp, value);
			}
		}
	}
	fclose(fr2);
	printf("Finish reading PMI file\n\n");

	HashNode* pNode = NULL;
	Model* model = NULL;

	int i1;
	for(i1 = 0; i1 < numM; ++i1){
		sprintf(s_tmp, "%s/%s", dir, model_dir[i1]);
		model = read_model(s_tmp);
		gsl_vector* vector = gsl_vector_calloc(model->v);
		double **pmi = (double**)malloc(sizeof(double*) * model->I);
		for(i = 0; i < model->I; ++i){
			pmi[i] = (double*)malloc(sizeof(double) * (model->K_I[i] + 1));
		}
		for(d = 0; d < model->I; ++d){//printf("calculate %d group\n", d);
			pmi[d][model->K_I[d]] = 0;
			for (k = 0; k < model->K_I[0]; k++){
				pmi[d][k] = 0;
				gsl_matrix_get_row(vector, model->phi_i[d], k);
				sort(vector);
				for(i = 0;i < topk; i++){
					for (j = i + 1;j < topk; j++){
						sprintf(s_tmp, "%d,%d", (int)vget(vector, i), (int)vget(vector, j));
						pNode = hash_table_lookup(s_tmp);
						if(pNode){
							value = pNode->nValue;
						}else{
							value = 0;
						}
						pmi[d][k] += value;
						pmi[d][model->K_I[d]] += value;
					}
				}
				pmi[d][k] /= ((topk - 1) * topk / 2);
			}
			pmi[d][model->K_I[d]] /= (model->K_I[d] * (topk - 1) * topk / 2);
		}

		value = 0;
		totV = 0;
		for(i = 0; i < model->I; ++i){
			value += pmi[i][model->K_I[i]] * model->W_I[i];
			totV += model->W_I[i];
		}
		value /= totV;
		printf("PMIs for model %s:\n", model_dir[i1]);
		for(d = 0; d < model->I; ++d){
			printf("%d:", d);
			for(j = 0; j < model->K_I[d]; ++j){
				printf(" %lf", pmi[d][j]);
			}
			printf("==>ave: %lf\n", pmi[d][model->K_I[d]]);
		}
		printf("ave: %lf\n\n", value);
		for(i = 0; i < model->I; ++i){
			free(pmi[i]);
		}
		free(pmi);
		free_model(model);
		gsl_vector_free(vector);
	}
	hash_table_release();
}

void deleteList(vlist* list)
{
	vlist* ll = list;
	while(ll){
		list = list->next;
		free(ll);
		ll = list;
	}
}
vlist* addListNode(vlist* list, int v)
{
	vlist* l1 = list;
	if(!list){
		l1 = (vlist*)malloc(sizeof(vlist));
		l1->len = 1;
		l1->v = v;
		l1->next = NULL;
	}else{
		while(list->next){
			(list->len)++;
			list = list->next;
		}
		(list->len)++;
		list->next = (vlist*)malloc(sizeof(vlist));
		list->next->len = list->len;
		list->next->v = v;
		list->next->next = NULL;
	}
	return l1;
}
vlist* removeListNode(vlist* list, int i)
{
	vlist* h;
	assert(list != NULL);
	assert(list->len > i);
	if(i == 0){
		h = list;
		list = list->next;
		free(h);
	}else{
		int ii;
		vlist* h1 = list;
		for(ii = 1; ii < i; ii++){
			h1 = h1->next;
		}
		h = h1->next;
		h1->next = h1->next->next;
		free(h);
	}
	h = list;
	while(h){
		(h->len)--;
		h = h->next;
	}
	return list;
}
int getNodeV(vlist* list, int i)
{
	int ii;
	assert(list);
	assert(i < list->len);
	for(ii = 0; ii < i; ii++){
		list = list->next;
	}
	return list->v;
}

/*
 * Note: free estimator must before free corpus!!!
 */
void free_estimator(Estimator* est, Corpus* c, int do_mu, int do_phi)
{
	int i;

	if(do_mu){
		for(i = 0; i < c->I; ++i)
			gsl_matrix_free(est->mu[i]);
		free(est->mu);
	}

	if(do_phi){
		gsl_matrix_free(est->phi);
		for(i = 0; i < c->I; ++i)
			gsl_matrix_free(est->phi_i[i]);
		free(est->phi_i);
	}
	free(est);
}

/*
 * Note: free topic assignment must before free corpus!!!
 */
void free_assignment(Assignment* ass, Corpus* c)
{
	int i, j;

	for (i = 0; i < c->I; i++){
		for (j = 0; j < c->ndocs[i]; j++){
			free(ass->topic_ass[i][j]);
		}
		free(ass->topic_ass[i]);
	}
	free(ass->topic_ass);
	free(ass);
}
/*
 * Note free plda_ct must before free corpus !!!
 */
void free_cts(Cts* cts, Corpus* c, int k)
{
	int i, j;

	free(cts->B);
	free(cts->sum_T);
	for(i = 0; i < k; ++i){
		free(cts->T_KQ[i]);
	}
	free(cts->T_KQ);
	for(i = 0; i < c->I; ++i){
		free(cts->N[i]);
		free(cts->top_cts[i].M);
		free(cts->top_cts[i].T_K);
		free(cts->top_cts[i].T_W);
		for(j = 0; j < c->ndocs[i]; ++j){
			free(cts->n[i][j]);
		}
		for(j = 0; j < k; ++j){
			free(cts->top_cts[i].m[j]);
			free(cts->top_cts[i].t[j]);
		}
	}
	free(cts->top_cts);
	free(cts->N);
	free(cts->n);
	free(cts);
}
void free_model(Model* model)
{
	int i, w;
	free(model->a);
	free(model->b);
	free(model->alpha_sum);
	gsl_vector_free(model->gamma);
	gsl_matrix_free(model->phi);
	if(model->rest_P){
		gsl_matrix_free(model->rest_P);
	}
	for(i = 0; i < model->I; ++i){
		gsl_matrix_free(model->phi_i[i]);

		if(model->P_k){
			gsl_vector_free(model->P_k[i]);
		}
		gsl_vector_free(model->alpha[i]);

		if(model->P_size){
			free(model->P_size[i]);
		}
		for(w = 0; w < model->W_I[i]; ++w){
			if(model->P){
				free(model->P[i][w]);
			}
			if(model->P_ind){
				free(model->P_ind[i][w]);
			}
		}
		if(model->P){
			free(model->P[i]);
		}
		if(model->P_ind){
			free(model->P_ind[i]);
		}
		int k;
		for(k = 0; k < model->K_I[0]; ++k){
			for(w = 0; w < model->W_I[i]; ++w){
				if(model->q){
					free(model->q[i][k][w]);
				}
				if(model->vL){
					deleteList(model->vL[i][k][w]);
				}
			}
			if(model->q){
				free(model->q[i][k]);
			}
			if(model->vL){
				free(model->vL[i][k]);
			}
		}
		if(model->q){
			free(model->q[i]);
		}
		if(model->vL){
			free(model->vL[i]);
		}
	}
	if(model->q){
		free(model->q);
	}
	if(model->vL){
		free(model->vL);
	}
	free(model->K_I);
	free(model->phi_i);
	free(model->W_I);
	if(model->P){
		free(model->P);
	}
	if(model->P_ind){
		free(model->P_ind);
	}
	if(model->P_size){
		free(model->P_size);
	}
	if(model->P_k){
		free(model->P_k);
	}
	free(model->alpha);

	if(model->P_size_all){
		free(model->P_size_all);
	}
	for(i = 0; i < model->v; ++i){
		if(model->P_all){
			free(model->P_all[i]);
		}
		if(model->P_idx_all){
			free(model->P_idx_all[i]);
		}
	}
	if(model->P_all){
		free(model->P_all);
	}
	if(model->P_idx_all){
		free(model->P_idx_all);
	}

	free(model);
}

/*
 *
 */
void free_all(Corpus* c, Model* model, Cts* cts, Assignment* ass, vocabulary* v)
{
	free_cts(cts, c, model->K_I[0]);
	free_model(model);
	free_assignment(ass, c);
	free_corpus(c);
	free_vocabulary(v);
	printf("--- All space is freed!!!\n");
}

