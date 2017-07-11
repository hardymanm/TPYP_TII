/*
 * est_gibbs.c
 *
 */

#include <string.h>
#include <math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include "gibbs.h"

extern double* vector;

static int sample_mul(double* q, int len)
{
	double ss, sum = 0;
	int i;
	for(i = 0; i < len; i++){
		sum += q[i];
	}
	sum *= next_uniform();
	ss = 0;
	for(i = 0; i < len; i++){
		ss += q[i];
		if(ss > sum && (q[i] != 0)){
			return i;
		}
	}
	return -1;
}

int remove_old_topic(int i, int d, int w, int k, Cts* cts, Model* model)
{
	int current_cust = cts->top_cts[i].m[k][w];
	int current_t = cts->top_cts[i].t[k][w];
	int v, r;
	if((double)current_t / (double)current_cust > next_uniform()){
		r = 1;
	}else{
		r = 0;
	}

	if(current_t == 1 && current_cust > 1){
		if(r == 1){
			return 0;
		}else{
			cts->top_cts[i].m[k][w]--;
			cts->top_cts[i].M[k]--;
			cts->n[i][d][k]--;
			cts->N[i][d]--;
		}
	}else if(current_t == 1 && current_cust == 1){
		if(r == 0){
			return 0;
		}else{
			int t;
			cts->top_cts[i].m[k][w]--;
			cts->top_cts[i].M[k]--;
			cts->n[i][d][k]--;
			cts->N[i][d]--;
			t = (int)(next_uniform() * cts->top_cts[i].t[k][w]);
			cts->top_cts[i].t[k][w]--;
			cts->top_cts[i].T_K[k]--;
			cts->top_cts[i].T_W[w]--;
			cts->sum_T[k]--;
			v = getNodeV(model->vL[i][k][w], t);
			model->q[i][k][w][v]--;
			cts->T_KQ[k][model->P_ind[i][w][v]]--;
			model->vL[i][k][w] = removeListNode(model->vL[i][k][w], t);
		}
	}else{
		cts->top_cts[i].m[k][w]--;
		cts->top_cts[i].M[k]--;
		cts->n[i][d][k]--;
		cts->N[i][d]--;
		if(r == 1){
			int t;
			t = (int)(next_uniform() * cts->top_cts[i].t[k][w]);
			cts->top_cts[i].t[k][w]--;
			cts->top_cts[i].T_K[k]--;
			cts->top_cts[i].T_W[w]--;
			cts->sum_T[k]--;
			v = getNodeV(model->vL[i][k][w], t);
			model->q[i][k][w][v]--;
			cts->T_KQ[k][model->P_ind[i][w][v]]--;
			model->vL[i][k][w] = removeListNode(model->vL[i][k][w], t);
		}
	}
	if(cts->top_cts[i].t[k][w] == 0){
		model->vL[i][k][w] = NULL;
	}
	return 1;
}

void add_new_topic (int i, int d, int w, int k, int v, Cts* cts, Model* model)
{assert(cts->top_cts[i].m[k][w] >= cts->top_cts[i].t[k][w]);
	cts->top_cts[i].m[k][w]++;
	cts->top_cts[i].M[k]++;
	cts->n[i][d][k]++;
	cts->N[i][d]++;

	v--;
	if(v >= 0){
		model->vL[i][k][w] = addListNode(model->vL[i][k][w], v);
		cts->top_cts[i].t[k][w]++;
		cts->top_cts[i].T_K[k]++;
		cts->top_cts[i].T_W[w]++;
		cts->sum_T[k]++;
		model->q[i][k][w][v]++;
		cts->T_KQ[k][model->P_ind[i][w][v]]++;
	}
}

static int sample_topic(int i, int d, int l, int w, Assignment* ass, Cts* cts, Model* model,
		Corpus* c)
{
	int old_topic, new_topic = -1;;
	int v, k, KK, sucess, idx;
	KK = model->K_I[i];

	old_topic = ass->topic_ass[i][d][l];
	sucess = remove_old_topic(i, d, w, old_topic, cts, model);
	if(sucess == 1){
		idx = 0;
		for(k = 0; k < KK; ++k){
			if(cts->top_cts[i].m[k][w] == 0){
				vector[idx] = 0;
			}else{
				vector[idx] = (vget(model->alpha[i], k) + (double)cts->n[i][d][k]) / (cts->B[k] + cts->top_cts[i].M[k]);
				vector[idx] *= (((double)(cts->top_cts[i].m[k][w] - cts->top_cts[i].t[k][w]) + 1) / (cts->top_cts[i].m[k][w] + 1));
				vector[idx] *= exp(stirling(cts->top_cts[i].m[k][w] + 1, cts->top_cts[i].t[k][w], model->a[k], i)
						- stirling(cts->top_cts[i].m[k][w], cts->top_cts[i].t[k][w], model->a[k], i));
			}
			idx++;
			for(v = 0; v < model->maxV; v++){
				if(v >= model->P_size[i][w]){
					vector[idx] = 0;
				}else{
					vector[idx] = model->P[i][w][v] * (vget(model->alpha[i], k) + cts->n[i][d][k]);
					vector[idx] *= (cts->B[k] + model->a[k] * cts->top_cts[i].T_K[k]) / (cts->B[k] + cts->top_cts[i].M[k]);
					vector[idx] *= ((double)(cts->top_cts[i].t[k][w] + 1)) / (cts->top_cts[i].m[k][w] + 1);
					vector[idx] *= (vget(model->gamma, model->P_ind[i][w][v]) + cts->T_KQ[k][model->P_ind[i][w][v]])
							/ (model->gamma_sum + cts->sum_T[k]);
					vector[idx] *= exp(stirling(cts->top_cts[i].m[k][w] + 1, cts->top_cts[i].t[k][w] + 1, model->a[k], i)
							- stirling(cts->top_cts[i].m[k][w], cts->top_cts[i].t[k][w], model->a[k], i));
				}
				idx++;
			}
		}
		idx = sample_mul(vector, idx);
		new_topic = idx / (model->maxV + 1);
		k = idx % (model->maxV + 1);
		add_new_topic (i, d, w, new_topic, k, cts, model);
		ass->topic_ass[i][d][l] = new_topic;
	}
	return new_topic;
}

void spyp_e_gibbs(Model* model, Cts* cts, Assignment* ass, Corpus* c)
{
	int i, l, d;

	for(i = 0; i < c->I; i++){
		for(d = 0; d < c->ndocs[i]; ++d){
			for(l = 0; l < c->docs[i][d].paras[0].total; ++l){
				sample_topic(i, d, l, c->docs[i][d].paras[0].words_local[l], ass,
						cts, model, c);
			}
		}

	}//printf("m1 = %d, t1 = %d, m2 = %d, t2 = %d\n", cts->top_cts[0].m[0][0], cts->top_cts[0].t[0][0], cts->top_cts[0].m[0][1], cts->top_cts[0].t[0][1]);
}
