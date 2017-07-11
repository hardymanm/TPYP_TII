/*
 * seqTM.c
 *
 */

#include <time.h>
#include <sys/types.h>
#include <sys/stat.h>

#include <gsl/gsl_vector.h>
#include <gsl/gsl_sf.h>
#include <gsl/gsl_math.h>
#include <unistd.h>

#include "seqTM.h"
#include "model.h"
#include "gibbs.h"
#include "util.h"
#include "params.h"

#define MAX_M 1000
#define SAMPLES 30
extern dpyp_params PARAMS;

double* vector;

double* lglik;

/*
 * Compute mu and nu
 *
 */
static void compute_mean_estimator(Estimator* est1, Estimator* est2, Model* model,
		Corpus* c, int n)
{
	int i, d, k;
	double mean;
	double val;

	for(i = 0; i < model->I; ++i){
		for(d = 0; d < c->ndocs[i]; d++){
			for(k = 0; k < model->K_I[i]; k++){
				mean = mget(est1->mu[i], d, k);
				val = mget(est2->mu[i], d, k);
				mean = (mean*(n -1) + val)/n;
				mset(est1->mu[i], d, k, mean);
			}
		}
	}
}

/*
 * Train the model
 *
 */
static void est_gibbs(Model* model, Cts* cts,
		Assignment* ass, Corpus* c, Corpus* c_te, vocabulary* v, char* root, int pred)
{
	int i, ite;
	int samples;
	int begin_k;
	clock_t t1, t2;
	double lg_lihood;
	char str[BUFSIZ];
	Corpus* c1 = NULL;
	Corpus* c2 = NULL;
	Assignment* ass_te = NULL;
	Cts* cts_te = NULL;

	Estimator*  est1 = new_estimator(c, model->I, model->K_I, model->W_I, model->v, 1, 1); // save mean value
	Estimator*  est2 = new_estimator(c, model->I, model->K_I, model->W_I, model->v, 1, 0);
	Estimator*  est_te = NULL;
	ite = 1;
	samples = 0;

	/*********** test data *************/
	if(c_te){
		if(pred == 0){
			c1 = (Corpus*)malloc(sizeof(Corpus));
			c2 = (Corpus*)malloc(sizeof(Corpus));
			split_corpus(c_te, c1, c2);
		}else{
			c1 = c_te;
		}
		est_te = new_estimator(c1, model->I, model->K_I, model->W_I, model->v, 1, 0);
		ass_te = new_assignment(c1);
		cts_te = new_cts(model->I, model->K_I, model->W_I, model->v, c1);
		init_state_gibbs(model, cts_te, ass_te, c1, 1, 0);
	}
	/*********************************/

	do{
		t1 = clock();
		spyp_e_gibbs(model, cts, ass, c);
		t2 = clock() - t1;
		t1 = clock();
		/************ test data **********/
		if(c_te){
			if(pred == 0){
				held_out_gibbs(model, cts_te, ass_te, c1);
			}
		}
		/********************************/
		t2 += (clock() - t1);
		printf(">>> running time for %d-iteration: %lf seconds\n",  ite, (double)t2/CLOCKS_PER_SEC);

		//printf("++++++++++++++++++++ SPYP EST model log likelihood for iteration %d = %lf\n", ite, lg_lihood);

		if(ite > PARAMS.burn_in && ite % PARAMS.sampling_lag == 0){
			samples++;
			//sprintf(str, "%s/SPYP_%d", root, samples);
			//mkdir(str, S_IRUSR | S_IWUSR | S_IXUSR);
			compute_phi(cts, model, c);
			//save_model(model, c, cts, model->phi, model->phi_i, str);

			compute_topic_proportion(est2, cts, model, c);
			compute_mean_estimator(est1, est2, model, c, samples);
			compute_mean_phi(est1->phi, est1->phi_i, model, samples);

			//printf("****** save top words and trained model ******\n");
			//begin_k = 0;
			/****************
			 * Global top words
			 */
			//sprintf(str, "%s/top_%d_words_%d.txt", root, PARAMS.top_words[0], ite);
			//print_top_words(PARAMS.top_words[0], begin_k, begin_k + model->K_I[0] - 1, model->phi, v, str);
			/******************
			 * Local top words
			 */
			//sprintf(str, "%s/top_%d_words_%d", root, PARAMS.top_words[0], ite);
			//print_local_top_words(PARAMS.top_words, begin_k, begin_k + model->K_I[0] - 1, model->phi_i,
			//		c, v, str);

			//sprintf(str, "word_asso_%d.txt", samples);
			//print_word_asso(model, v, str);
		}

		if(ite % 3 == 0){

			lg_lihood = spyp_marg_lglihood1(c, cts, model, est2->mu, 1);
			printf("\n++++++++++++ SPYP train log-likelihood +++++++++++++++\n");
			printf("group-1: %lf", lglik[0]);
			for(i = 1; i < model->I; i++){
				printf(", group-%d: %lf", i + 1, lglik[i]);
			}
			printf(" ===> ave: %lf\n\n", lg_lihood);

			if(c_te){
				if(pred == 0){
					lg_lihood = spyp_marg_lglihood1(c2, cts_te, model, est_te->mu, 0);
					printf("\n++++++++++++ SPYP test perplexity +++++++++++++++\n");
					printf("group-1: %lf", gsl_sf_exp(-lglik[0] / c2->total[0]));
					lg_lihood = lglik[0];
					for(i = 1; i < model->I; i++){
						printf(", group-%d: %lf", i + 1, gsl_sf_exp(-lglik[i] / c2->total[i]));
						lg_lihood += lglik[i];
					}
					printf(" ===> ave: %lf\n\n", gsl_sf_exp(-lg_lihood / c2->totalw));
				}else if(ite > PARAMS.burn_in){
					lg_lihood = left2right_predict(model, cts_te, ass_te, c1, 4, 2);
					printf("left to right prediction accuracy: %lf\n\n", lg_lihood);
				}
			}

			//alpha_opt(cts, model, c->ndocs, model->K_I[0]);
			//alpha_opt_newton(cts, model, c->ndocs, model->K_I[0]);
			//sample_alpha(cts, model, c);
			printf("opt_alpha is done!!!\nalpha = ");
			for(i = 0; i < model->I; ++i){
				printf("%f ", gsl_vector_get(model->alpha[i], 0));
			}
			printf("\n");
			//gamma_opt(cts, model, c);
			//gamma_opt_newton(cts, model, c);
			printf("gamma = %lf\n", vget(model->gamma, 0));

			if(ite > 0){//50){
#ifdef SAB
				sample_ab(2, cts, model, c);
				update_stirling_table(model->a, model->I);
#else
				sample_b(2, cts, model, c);
				//sample_b_all(2, cts, model, c);
#endif
			}
		}

		ite++;
	}while(ite <= PARAMS.gibbs_max_iter);

	sprintf(str, "%s/SPYP_final", root);
	mkdir(str, S_IRUSR | S_IWUSR | S_IXUSR);
	//save_model(model, c, cts, est1->phi, est1->phi_i, str);
	save_model(model, c, cts, est1->phi, model->phi_i, str);

	save_customer_table(cts, model->I, model->K_I[0], model->v, str);

	save_estimator(est1, c, str, 1, 0);
	sprintf(str, "%s/SPYP_est_final", root);
	mkdir(str, S_IRUSR | S_IWUSR | S_IXUSR);
	printf_estimator(est1, c, str, 1, 0);

	begin_k = 0;
	/****************
	 * Global top words
	 */
	sprintf(str, "%s/top_%d_words_final.txt", root, PARAMS.top_words[0]);
	print_top_words(PARAMS.top_words[0], begin_k, begin_k + model->K_I[0] - 1, model->phi, cts, v, str);
	/*******************
	 * Local top words
	 */
	sprintf(str, "%s/top_%d_words_final", root, PARAMS.top_words[0]);
	print_local_top_words(PARAMS.top_words, begin_k, begin_k + model->K_I[0] - 1, model->phi_i,
			cts, c, v, str);
	save_topic_assignmnet(c, ass, root);

	free_estimator(est1, c, 1, 1);
	free_estimator(est2, c, 1, 0);

	if(c_te){
		free_estimator(est_te, c1, 1, 0);
		free_assignment(ass_te, c1);
		free_cts(cts_te, c1, model->K_I[0]);
		if(pred == 0){
			free_corpus(c1);
			free_corpus(c2);
		}
	}
}

void estimate(Corpus* c, Corpus* c_te, vocabulary* v, char* root, char* dir_p, int init, char* z_root)
{
	Model* model;
	Cts* cts;
	Assignment* ass;

	initial_rng();

#ifdef IP
	dir_p = NULL;
#endif

	model = random_init(PARAMS.I, PARAMS.K_i, c->W_i, v->size, PARAMS.alpha, PARAMS.gamma,
			PARAMS.a, PARAMS.b, c->l2g_vocabulary, dir_p);

	cts = new_cts(model->I, model->K_I, model->W_I, model->v, c);
	if(init){
		ass = read_topic_assignmnet(c, z_root);
	}else{
		ass = new_assignment(c);
	}
	/*
	 * pre-compute tables:
	 * 		1. stirling table
	 */
	if(c->W_i[0] / 2 < MAX_M){
		make_stirling_table(c->W_i[0] / 2, c->W_i[0] / 4, model->a, model->I);
	}else{
		make_stirling_table(c->W_i[0] / 2, MAX_M, model->a, model->I);
	}

	printf(">>>>>> Begin Gibbs sampling iteration for SPYP Model ......\n");
	if(init)
		init_state_gibbs(model, cts, ass, c, 0, 1);
	else
		init_state_gibbs(model, cts, ass, c, 0, 0);

	int i, w;
	model->maxV = 0;
	for(i = 0; i < model->I; i++){
		for(w = 0; w < model->W_I[i]; w++){
			if(model->P_size[i][w] > model->maxV){
				model->maxV = model->P_size[i][w];
			}
		}
	}
	vector = (double*)malloc(sizeof(double) * (model->maxV + 1) * model->K_I[0]);

	lglik = (double*)malloc(sizeof(double) * model->I);

	est_gibbs(model, cts, ass, c, c_te, v, root, 0);
	printf(">>>>>> End Gibbs sampling iteration \n");

	free(lglik);
	free_stirling_table(model->I);
	free_rng();
	free_all(c, model, cts, ass, v);
	free(vector);
	if(c_te){
		free_corpus(c_te);
	}
}

/*
 *  Model log likelihood of the sequenctial model
 *
 */
double spyp_model_lglihood(Corpus* c, Cts* cts, Model* model)
{
	int i, d, k, w, v;
	double lglihood = 0;
	double* alpha_lg, alpha_sum_lg, gamma_sum_lg;

	for(i = 0; i < model->I; ++i){
		for(k = 0; k < model->K_I[i]; ++k){
			lglihood += log_poch_sym(cts->B[k], model->a[k], cts->top_cts[i].T_K[k])
					- log_poch_sym(model->b[k], 1, cts->top_cts[i].M[k]);
			for(w = 0; w < model->W_I[i]; ++w){
				lglihood += stirling(cts->top_cts[i].m[k][w], cts->top_cts[i].t[k][w], model->a[k], i);
				if(!gsl_finite(lglihood)){
					printf("m = %d, k = %d\n", cts->top_cts[i].m[k][w], cts->top_cts[i].t[k][w]);
					exit(0);
				}
				for(v = 1; v <= cts->top_cts[i].t[k][w]; ++v){
					lglihood += gsl_sf_log(v);
				}
				for(v = 1; v <= cts->top_cts[i].m[k][w] - cts->top_cts[i].t[k][w]; ++v){
					lglihood += gsl_sf_log(v);
				}
				for(v = 1; v <= cts->top_cts[i].m[k][w]; ++v){
					lglihood -= gsl_sf_log(v);
				}
			}
		}
		for(w = 0; w < model->W_I[i]; w++){
			for(v = 0; v < model->P_size[i][w]; v++){
				d = 0;
				for(k = 0; k < model->K_I[i]; k++){
					d += model->q[i][k][w][v];
				}
				lglihood += d * gsl_sf_log(model->P[i][w][v]);
			}
		}
		alpha_lg = (double*)malloc(sizeof(double) * model->K_I[i]);
		for(k = 0; k < model->K_I[i]; ++k)
			alpha_lg[k] = gsl_sf_lngamma(vget(model->alpha[i], k));
		alpha_sum_lg = gsl_sf_lngamma(model->alpha_sum[i]);
		for(d = 0; d < c->ndocs[i]; ++d){
			lglihood += alpha_sum_lg - gsl_sf_lngamma(model->alpha_sum[i] + cts->N[i][d]);
			for(k = 0; k < model->K_I[i]; ++k){
				lglihood += gsl_sf_lngamma(vget(model->alpha[i], k) + cts->n[i][d][k]) - alpha_lg[k];
			}
		}
		free(alpha_lg);
	}
	gamma_sum_lg = gsl_sf_lngamma(model->gamma_sum);
	for(k = 0; k < model->K_I[0]; ++k){
		lglihood += gamma_sum_lg - gsl_sf_lngamma(model->gamma_sum + cts->sum_T[k]);
		for(v = 0; v < model->v; ++v){
			lglihood += gsl_sf_lngamma(vget(model->gamma, v) + cts->T_KQ[k][v]) - gsl_sf_lngamma(vget(model->gamma, v));
		}
	}
	if(!gsl_finite(lglihood)){
		fprintf(stderr, "Error: DPYP model log likelihood is NaN or infinite!!!\n");
		exit(1);
	}
	return lglihood;
}

double spyp_marg_lglihood(Corpus* c, Cts* cts, Assignment* ass, Model* model, int inf)
{
	int i, k, d, l, v;
	double lglik_ave = 0;

	if(inf == 0){
		compute_phi(cts, model, c);
	}
	for(i = 0; i < model->I; i++){
		lglik[i] = 0;
		for(d = 0; d < c->ndocs[i]; ++d){
			for(l = 0; l < c->docs[i][d].paras[0].total; ++l){
				k = ass->topic_ass[i][d][l];
				if(inf == 1){
					v = c->docs[i][d].paras[0].words[l];
				}else{
					v = c->docs[i][d].paras[0].words_local[l];
				}
				lglik[i] += gsl_sf_log(mget(model->phi_i[i], k, v));
			}
		}
		lglik_ave += lglik[i] * c->total[i] / c->totalw;
	}
	return lglik_ave;
}

double spyp_marg_lglihood1(Corpus* c, Cts* cts, Model* model, gsl_matrix** mu, int train)
{
	int i, k, d, l, w;
	double lglik_ave = 0;

	if(train == 1){
		compute_phi(cts, model, c);
	}
	compute_mu(cts, model, c, mu);
	for(i = 0; i < c->I; i++){
		lglik[i] = 0;
		for(d = 0; d < c->ndocs[i]; d++){
			for(l = 0; l < c->docs[i][d].paras[0].total; l++){
				double val = 0;
				w = c->docs[i][d].paras[0].words_local[l];
				for(k = 0; k < model->K_I[i]; k++){
					val += fabs(mget(model->phi_i[i], k, w)) * mget(mu[i], d, k);
				}
				lglik[i] += gsl_sf_log(val);
			}
		}
		lglik_ave += lglik[i] * c->total[i] / c->totalw;
	}
	return lglik_ave;
}
