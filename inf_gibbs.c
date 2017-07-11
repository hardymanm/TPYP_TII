/*
 * inf_Gibbs.c
 *
 */

#include "gibbs.h"
#include "params.h"


extern dpyp_params PARAMS;
extern double* vector;

static int held_out_sample_topic(int i, int d, int l, int w, Assignment* ass, Cts* cts, Model* model)
{
	int old_topic, new_topic;
	int k, KK;
	double sum_v;

	KK = model->K_I[i];
	old_topic = ass->topic_ass[i][d][l];
	cts->n[i][d][old_topic] -= 1;

	for(k = 0; k < KK; ++k){
		double r = (vget(model->alpha[i], k) + cts->n[i][d][k]) * mget(model->phi_i[i], k, w);
		vector[k] = r;
	}

	sum_v = 0;
	for(k = 0; k < KK; ++k)
		sum_v += vector[k];
	for(k = 0; k < KK; ++k)
		vector[k] = vector[k] / sum_v;
	new_topic = next_discrete_normalised(vector, KK);

	ass->topic_ass[i][d][l] = new_topic;
	cts->n[i][d][new_topic] += 1;

	return new_topic;
}

static int held_out_sample_topic_pred(int i, int phi_i, int d, int l, int w, Assignment* ass, Cts* cts, Model* model)
{
	int old_topic, new_topic;
	int k, KK;
	double sum_v;

	KK = model->K_I[phi_i];
	old_topic = ass->topic_ass[i][d][l];
	cts->n[i][d][old_topic] -= 1;

	for(k = 0; k < KK; ++k){
		double r = (vget(model->alpha[phi_i], k) + cts->n[i][d][k]) * mget(model->phi_i[phi_i], k, w);
		vector[k] = r;
	}

	sum_v = 0;
	for(k = 0; k < KK; ++k)
		sum_v += vector[k];
	for(k = 0; k < KK; ++k)
		vector[k] = vector[k] / sum_v;
	new_topic = next_discrete_normalised(vector, KK);

	ass->topic_ass[i][d][l] = new_topic;
	cts->n[i][d][new_topic] += 1;

	return new_topic;
}

void held_out_gibbs(Model* model, Cts* cts, Assignment* ass, Corpus* c)
{
	int i, d, l;

	for(i = 0; i < c->I; i++){
		for(d = 0; d < c->ndocs[i]; ++d){
			for(l = 0; l < c->docs[i][d].paras[0].total; ++l){
				held_out_sample_topic(i, d, l, c->docs[i][d].paras[0].words_local[l], ass, cts, model);
			}
		}
	}
}

double left2right_predict(Model* model, Cts* cts, Assignment* ass, Corpus* c, int nSamples, int burnin)
{
	int i, ii, k, r, d, l, lp, pl, pl_m, idx, d_t;
	double correct;

	if(nSamples <= burnin){
		nSamples = burnin + 1;
	}

	d_t = 0;
	correct = 0;
	for(i = 0; i < c->I; i++){
		for(d = 0; d < c->ndocs[i]; ++d){
			d_t++;
			idx = 0;
			pl_m = 0;
			for(ii = 0; ii < model->I; ii++){

				pl = 0;
				for(lp = 0; lp < c->docs[i][d].paras[0].total; ++lp){
					for(k = 0; k < model->K_I[ii]; k++){
						cts->n[i][d][k] = 0;
					}
					for(l = 0; l < lp; l++){
						cts->n[i][d][ass->topic_ass[i][d][l]]++;
					}

					for(r = 0; r < nSamples; r++){
						for(l = 0; l < lp; l++){
							held_out_sample_topic_pred(i, ii, d, l, c->docs[i][d].paras[0].words_local[l], ass, cts, model);
						}
						if(r >= burnin){
							for(k = 0; k < model->K_I[ii]; k++){
								pl += log((vget(model->alpha[ii], k) + cts->n[i][d][k])
										/(model->alpha_sum[ii] + lp) * mget(model->phi_i[ii], k, c->docs[i][d].paras[0].words_local[lp]));
							}
						}
					}
				}
				if(ii == 0){
					pl_m = pl;
					idx = ii;
				}else{
					if(pl > pl_m){
						pl_m = pl;
						idx = ii;
					}
				}
			}
			if(idx == i){
				correct++;
			}
		}
	}
	return correct / d_t;
}
