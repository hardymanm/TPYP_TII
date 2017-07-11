/*
 * gibbs.c
 *
 */
#include "gibbs.h"

/*
 * uniformly initialize the topic assignment of each word
 */

void init_state_gibbs(Model* model, Cts* cts, Assignment* ass, Corpus* c, int est_inf, int init_z)
{
	int i, j, l, k, q, v, w;
	time_t seed;
	gsl_rng* r;

	printf(">>>>>> Begin initialize the states ......\n");
	cts->I = model->I;
	for(i = 0; i < model->K_I[0]; ++i)
		cts->B[i] = model->b[i];
	for(k = 0; k < model->K_I[0]; ++k)
		cts->sum_T[k] = 0;
	for(i = 0; i < c->I; ++i){
		for(j = 0; j < c->ndocs[i]; ++j){
			cts->N[i][j] = 0;
			for(k = 0; k < model->K_I[i]; ++k){
				cts->n[i][j][k] = 0;
			}
		}
		for(k = 0; k < model->K_I[i]; ++k){
			cts->top_cts[i].M[k] = 0;
			cts->top_cts[i].T_K[k] = 0;
			for(w = 0; w < model->W_I[i]; ++w){
				cts->top_cts[i].m[k][w] = 0;
				cts->top_cts[i].t[k][w] = 0;
			}
		}
		for(w = 0; w < model->W_I[i]; ++w){
			cts->top_cts[i].T_W[w] = 0;
		}
	}
	r = gsl_rng_alloc(gsl_rng_taus);
	time(&seed);
//	seed = 1115574245;
	gsl_rng_set(r, (long)seed);

	for (i = 0; i < c->I; i++){
		for (j = c->ndocs[i] - 1; j >= 0; j--){
			for (l = 0; l < c->docs[i][j].paras[0].total; l++){
				if(init_z){
					q = ass->topic_ass[i][j][l];
				}else{
					q = (int)floor(gsl_rng_uniform(r) * model->K_I[i]);
					ass->topic_ass[i][j][l] = q;
				}
				cts->n[i][j][q] += 1;
				cts->N[i][j] += 1;
				if(est_inf == 0){
					cts->top_cts[i].m[q][c->docs[i][j].paras[0].words_local[l]] += 1;
					cts->top_cts[i].M[q] += 1;
				}
			}
		}
		if(est_inf == 0){
			for(k = 0; k < model->K_I[i]; ++k){
				for(w = 0; w < model->W_I[i]; ++w){
					if(cts->top_cts[i].m[k][w]){
						cts->top_cts[i].t[k][w] = 1;
						cts->top_cts[i].T_K[k] += 1;
						cts->top_cts[i].T_W[w] += 1;
						cts->sum_T[k]++;
					}
				}
			}
		}
	}
	if(est_inf == 0){
		int v1;
		for(i = 0; i < c->I; i++){
			for(k = 0; k < model->K_I[i]; k++){
				for(w = 0; w < model->W_I[i]; w++){
					model->vL[i][k][w] = NULL;
					if(cts->top_cts[i].m[k][w] > 0){
						for(v = 0; v < model->P_size[i][w]; v++){
							if(c->l2g_vocabulary[i][w] == model->P_ind[i][w][v]){
								model->q[i][k][w][v] = cts->top_cts[i].t[k][w];
								v1 = v;
								break;
							}
						}
						assert(v < model->P_size[i][w]);
						model->vL[i][k][w] = addListNode(model->vL[i][k][w], v1);
					}
				}
			}
		}
		for(k = 0; k < model->K_I[0]; ++k){
			for(v = 0; v < model->v; ++v){
				cts->T_KQ[k][v] = 0;
			}
		}
		for(k = 0; k < model->K_I[0]; ++k){

			for(i = 0; i < model->I; ++i){
				for(w = 0; w < model->W_I[i]; ++w){
					for(v = 0; v < model->P_size[i][w]; ++v){
						cts->T_KQ[k][model->P_ind[i][w][v]] += model->q[i][k][w][v];
					}
				}
			}
		}
	}

	gsl_rng_free(r);
	printf("The seed for init_gibbs: %ld\n", (long)seed);
	printf("<<<<<< End initialization \n");
}

/*
 * compute matrix phi
 */
void compute_phi(Cts* cts, Model* model, Corpus* c)
{
	int k, v, i, w;

	/********************
	 * Global topic-word matrix computation
	 */
	for(k = 0; k < model->K_I[0]; ++k){
		for(v = 0; v < model->phi->size2; ++v){
			mset(model->phi, k, v, (cts->T_KQ[k][v] + vget(model->gamma, v)) / (cts->sum_T[k] + model->gamma_sum));
		}
	}
	/*****************
	 * Local topic-word matrix computation
	 */
	for(i = 0; i < model->I; ++i){
		gsl_vector *phi_tmp = gsl_vector_calloc(model->phi_i[i]->size2);
		for(k = 0; k < model->K_I[i]; k++){
			gsl_vector_set_all(phi_tmp, 0);
			for(w = 0; w < phi_tmp->size; ++w){
				for(v = 0; v < model->P_size[i][w]; ++v){
					vinc(phi_tmp, w, mget(model->phi, k, model->P_ind[i][w][v]) * model->P[i][w][v]);
				}
			}

			for(w = 0; w < phi_tmp->size; ++w){
				mset(model->phi_i[i], k, w,
						((cts->top_cts[i].m[k][w] - model->a[k] * cts->top_cts[i].t[k][w])
								+ vget(phi_tmp, w) * (cts->top_cts[i].T_K[k] * model->a[k] + model->b[k]))
								//+ mget(model->phi, k, c->l2g_vocabulary[i][w]) * (cts->top_cts[i].T_K[k] * model->a[k] + model->b[k]))
								/ (model->b[k] + cts->top_cts[i].M[k]));
			}
		}
		gsl_vector_free(phi_tmp);
	}
}

void compute_mean_phi(gsl_matrix* mphi, gsl_matrix**mphi_i, Model* model, int n)
{
	int i, q, w;
	double val;
	assert(mphi->size1 == model->phi->size1 && mphi->size2 == model->phi->size2);
	for(q = 0; q < model->phi->size1; q++){
		for(w = 0; w < model->phi->size2; w++){
			val = (mget(mphi, q, w)*(n - 1) + mget(model->phi, q, w))/n;
			mset(mphi, q, w, val);
		}
	}
	for(i = 0; i < model->I; ++i){
		for(q = 0; q < model->phi_i[i]->size1; q++){
			for(w = 0; w < model->phi_i[i]->size2; w++){
				val = (mget(mphi_i[i], q, w)*(n - 1) + mget(model->phi_i[i], q, w))/n;
				mset(mphi_i[i], q, w, val);
			}
		}
	}
}

void compute_mu(Cts* cts, Model* model, Corpus* c, gsl_matrix** mu)
{
	int i, ii, q;
	for(ii = 0; ii < model->I; ++ii){
		for(i = 0; i < c->ndocs[ii]; i++){
			for(q = 0; q < model->K_I[ii]; q++){
				double val = (vget(model->alpha[ii], q) + cts->n[ii][i][q])
						/(model->alpha_sum[ii] + cts->N[ii][i]);
				mset(mu[ii], i, q, val);
			}
		}
	}
}

void compute_topic_proportion(Estimator* estimator, Cts* cts, Model* model, Corpus* c)
{
	compute_mu(cts, model, c, estimator->mu);
}
